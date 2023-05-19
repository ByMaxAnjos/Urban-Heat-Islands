#=======================================================================================
#
# Title:       High Spatio-temporal Resolution of Urban Heat Islands
# Author:      Dr. Max Anjos (maxanjos@campus.ul.pt)
# Description: Interpolation of air temperarure. More details on the approach are available at:
#             https://github.com/ByMaxAnjos/Urban-Heat-Islands
# Data: 2.05.2023


#' @param air_temperature_data (required). Air temperature timseries (.csv format) from different local stations that should include latitude, longitude columun
#' @param cws_data (conditionally required).
#' @param weather (required) Hourly meteorological data from the German Weather Service (DWD package).
#' @return Air temperature summary, csv. table and sf multipolylines and raster maps
#' @examples
#=======================================================================================

# Load theses packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, sf, data.table, raster, tmap, osmdata, gstat, automap, openair, recipes, timetk)
library(tidyverse)
library(data.table)
library(tmap)
library(osmdata)
library(sf)
library(raster)
library(gstat)   # The most popular R-Package for Kriging (imho)
library(automap) # Automatize some (or all) parts of the gstat-workflow
library(openair)
library(recipes)
library(timetk)
library(terra)
library(stars)
library(tidymodels)
#Define your path
setwd("/Users/test/Documents/CO2CityMap/Berlin/Components/Building")

#================================================================
#Load air temperature data and pre-processing
#================================================================

#CWS data
#air_cws_2018 <- fread("inputs/data/cws_berlin_ta_level_o1_2018_July.csv")
#air_cws_2019 <- fread("nputs/data/cws_berlin_ta_level_o1_2019_July.csv")

#Pre=processing data
air_cws_2018 <- air_cws_2018 %>%
  rename(date=time,
         longitude = lon,
         latitude = lat,
         airT = ta_int) %>%
  dplyr::select(date, longitude, latitude, airT)

#UCON data
air_UCON <- fread("inputs/data/airT_UCON_2015_2022.csv") %>%
  rename(longitude = Longitude,
         latitude = Latitude)

#================================================================
#Define the region of int/project system in UTM
#================================================================

#get ROI
city <- "Berlin"
mycrs <- "+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs"

study_area <- osmdata::getbb(city, format_out = "sf_polygon", limit = 1)$multipolygon %>%
  st_as_sf() %>% # convert to sf object
  st_transform(crs = mycrs)
study_area <- st_make_valid(study_area)
qtm(study_area) #plot map

# study_area <- osmdata::getbb(city, format_out = "sf_polygon", limit = 1) %>% #orthewise try this one
#   st_as_sf() %>% # convert to sf object
#   st_transform(crs = mycrs)
# study_area <- st_make_valid(study_area)
# qtm(study_area) #plot map

#================================================================
#Load LCZ map
#================================================================
lcz_map <- raster("LCZ/lcz_berlin2.tif")
lcz_map <- raster::projectRaster(lcz_map, crs = mycrs)

#Building
build <- raster::raster("building/inputs/raster/build100.tif")
names(build) <- paste0("volume")


#define the period
#month <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
imonth <- c("feb")
iyear <- c(2019)
myday <- 2
idates <- expand.grid(imonth, iyear)

#================================================================
#Deploy the function for interpolation of air temperature
#================================================================
'
idates is the dataframe with the defined period.
air_df is the pre-processed air temperature data
spRes is the horizontal resolution in meters
impute if TRUE impute the missing air "NaN" crws data.
varmodel is the variogram model type. It can be "Exp", "Sph", "Gau" or "Mat".
isave if TRUE creates a new folder "ZCCM_Output" in you path and saves the rasters
'

AirLCZidw <- function(idates,
                           air_df = air_UCON,
                           roi = study_area,
                           spRes = 500,
                           lcz = lcz_map,
                           impute = FALSE,
                           IDW = TRUE,
                           evaluation = FALSE,
                           isave = FALSE) {
  imonth <- idates[1]
  iyear <- idates[2]

  if(impute ==TRUE) {
    air_df[air_df$airT == "NAN"] <- NA
    air_recipe <- recipe(date ~., data = air_df) %>%
      step_ts_impute(all_predictors(), -latitude, -longitude)
    air_model <- air_recipe %>%
      prep(air_df) %>%
      bake(air_df)
  } 

  #Pre=processing data
  if(is.null(myday)) {
    air_model <- air_df %>%
      openair::selectByDate(year = iyear, month = imonth) %>%
      na.omit()
  } else {
    air_model <- air_df %>%
      openair::selectByDate(year = iyear, month = imonth, day = myday) %>%
      na.omit()
  }
  
  # Convert to spatial pixel
  grid_stations <- raster::raster(raster::extent(roi), res = spRes)
  raster::values(grid_stations) <- 1:raster::ncell(grid_stations)
  raster::crs(grid_stations) <- mycrs

  iLCZ <-
    raster::crop(lcz, extent(roi)) %>%
    raster::mask(roi) %>%
    resample(grid_stations) #Resample
  names(iLCZ) <- "lcz"
  #Convert lcz_map to polygon
  lcz_shp <- terra::as.polygons(rast(iLCZ)) %>%
    st_as_sf() %>%
    st_transform(crs = mycrs)
  #Convert to list of polygon
  poly_list <- lapply(st_geometry(lcz_shp), function(x)
    st_sfc(x, crs = mycrs))
  #Calculate areas of each polygon
  lcz_areas <- lapply(1:length(poly_list), function(i)
    st_area(poly_list[[i]]))
  #Sample the number of points according to area of the polygon and convert to data.frame
  lcz_poi <- lapply(1:length(poly_list), function(i)
    st_sample(poly_list[[i]], size=100, prob=lcz_areas[[i]], method = "random", exact = FALSE) %>%
      as.data.frame())
  #Merge all dataframes and convert to sf points
  lcz_poi <- do.call(rbind.data.frame, lcz_poi) %>%
    st_as_sf() %>%
    st_transform(crs = mycrs)
  #Intersect lcz poi with lcz shp
  lcz_poi_get <- sf::st_intersection(lcz_poi, lcz_shp)

  #Downscale to hour
  iday <- air_model %>%
    mutate(iday = lubridate::day(date)) %>%
    distinct(iday, .keep_all = FALSE) %>%
    expand.grid()


  if(IDW == TRUE) {

    model_day <- function(iday) {

    myday <- iday[1]

    modelday <- air_model %>%
      mutate(day = lubridate::day(date)) %>%
      openair::selectByDate(day = myday)

    #Downscale to hour
    ihour <- modelday %>%
      mutate(ihour = hour(date)) %>%
      distinct(ihour, .keep_all = FALSE) %>%
      expand.grid()

    model_hour <- function(ihour) {

      myhour <- ihour[1]

      data_model <- modelday %>%
        mutate(hour = lubridate::hour(date)) %>%
        openair::selectByDate(hour = myhour)

        #Get station points
        shp_stations <- data_model %>%
          distinct(latitude, longitude, .keep_all = T) %>%
          dplyr::select(station, latitude, longitude, airT) %>%
          sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
          st_transform(crs = mycrs) 

        #Intersect shp stations with lcz shp
        lcz_stations <- sf::st_intersection(shp_stations, lcz_shp)
        #Create train and test sets
        lcz_stations_split <- initial_split(lcz_stations)
        lcz_stations_train <- training(lcz_stations_split)
        lcz_stations_test <- testing(lcz_stations_split)

        #Merge table to create a sf to model
       lcz_poi_mod <- inner_join(lcz_poi_get,lcz_stations_train %>% as_tibble() %>%
                                  dplyr::select(lcz, station, airT), by=c("lcz")) %>%
          group_by(station) %>%
          mutate(FID_station=cur_group_id()) %>%
          dplyr::select(-station, -lcz)

        # Convert LCZ map to starts
        lcz_stars <- st_as_stars(iLCZ, dimensions = "XY")
        train_mod = st_join(lcz_poi_mod, st_as_sf(lcz_stars)) %>%
          mutate(lcz = as.integer(lcz))
        train_mod = train_mod[!is.na(train_mod$lcz),]
        st_crs(lcz_stars) <- st_crs(train_mod)

        #IDW
          idw_mod = gstat(formula = airT ~ 1, data = train_mod)
          idw_map = predict(idw_mod, lcz_stars)
          idw_map = idw_map["var1.pred",,]
          mydate <- data_model %>% distinct(date, .keep_all = FALSE) %>% as.data.frame()
          mydate <- gsub("[: -]", "" , mydate$date, perl=TRUE)
          names(idw_map) <- paste0("IDW_airT_",mydate)
          idw_map <- raster(rast(idw_map))
          #return(idw_map)

          # Resample air map
          air_resample = raster::resample(idw_map, build)
          #raster::writeRaster(air_resample, paste0("Building/outputs/maps/UHI/", mydate, "_UHI.TIF"), format="GTiff", overwrite = TRUE)
          #Merge building with air raster
          air_build <- as_tibble(rasterToPoints(air_resample)) %>% set_names(c("x", "y", "airT"))
          build_df <- as_tibble(rasterToPoints(build))
          build_model <- inner_join(build_df, air_build, by= c("x", "y")) %>%
            mutate(hour = paste0(myhour))

          #Calculate CO2 emissions
          ECO2B <- ECO2build(build_model)
          ECO2T_ras <- ECO2B %>% dplyr::select(x, y, ECO2_micro)
          ECO2T_ras = raster::rasterFromXYZ(xyz = ECO2T_ras,crs = mycrs)
          names(ECO2T_ras) <- paste0("CO2B_", mydate)
          return(ECO2T_ras)
      #Salve the map

          if(isave == TRUE) {

            # Create a folder name using paste0
            folder <- paste0("ZCCM_Output/")

            # Check if the folder exists
            if (!dir.exists(folder)) {
              # Create the folder if it does not exist
              dir.create(folder)
            }

            #Salve the map
            file <- paste0(folder,iyear,imonth,myday,myhour,"residual_map.png")
            ggsave(iplot,file)

            file <- paste0(folder,iyear,imonth,myday,myhour,"AirT_map.png")
            mapview::mapshot(iplot2,file)

          }

         
        }

    MapHour <- apply(ihour, 1, model_hour)
    return(MapHour)

    }

    MapDay <- pbapply::pbapply(iday, 1, model_day)
    #MapDay <- unlist(MapDay)
    return(MapDay)

  }

  if(evaluation == TRUE) {

    model_day <- function(iday) {

      myday <- iday[1]

      modelday <- air_model %>%
        mutate(day = lubridate::day(date)) %>%
        openair::selectByDate(day = myday)

      #Downscale to hour
      ihour <- modelday %>%
        mutate(ihour = hour(date)) %>%
        distinct(ihour, .keep_all = FALSE) %>%
        expand.grid()

      model_hour <- function(ihour) {

        myhour <- ihour[1]

        data_model <- modelday %>%
          mutate(hour = lubridate::hour(date)) %>%
          openair::selectByDate(hour = myhour)

        #Get station points
        shp_stations <- data_model %>%
          distinct(latitude, longitude, .keep_all = T) %>%
          dplyr::select(station, latitude, longitude, airT) %>%
          sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
          st_transform(crs = mycrs) %>%
          sf::st_intersection(roi)

        #Intersect shp stations with lcz shp
        lcz_stations <- sf::st_intersection(shp_stations, lcz_shp)
        #Create train and test sets
        lcz_stations_split <- initial_split(lcz_stations)
        lcz_stations_train <- training(lcz_stations_split)
        lcz_stations_test <- testing(lcz_stations_split)

        #Merge table to create a sf to model
        lcz_poi_mod <- inner_join(lcz_poi_get,lcz_stations_train %>% as_tibble() %>%
                                    dplyr::select(lcz, station, airT), by=c("lcz")) %>%
          group_by(station) %>%
          mutate(FID_station=cur_group_id()) %>%
          dplyr::select(-station, -lcz)

        # Convert LCZ map to starts
        lcz_stars <- st_as_stars(iLCZ, dimensions = "XY")
        train_mod = st_join(lcz_poi_mod, st_as_sf(lcz_stars)) %>%
          mutate(lcz = as.integer(lcz))
        train_mod = train_mod[!is.na(train_mod$lcz),]
        st_crs(lcz_stars) <- st_crs(train_mod)

        #IDW
        idw_mod = gstat(formula = airT ~ 1, data = train_mod)
        idw_map = predict(idw_mod, lcz_stars)
        idw_map = idw_map["var1.pred",,]
        # mydate <- data_model %>% distinct(date, .keep_all = FALSE) %>% as.data.frame()
        # mydate <- gsub("[: -]", "" , mydate$date, perl=TRUE)
        names(idw_map) <- paste0("IDW_airT")

        #Evaluation
        idw_cv = gstat.cv(idw_mod, nfold = 5, verbose = FALSE)
        idw_cv = st_as_sf(idw_cv)

        #This is the RMSE value for the IDW interpolation with sampled points
        rmse_sampled= data_model %>% distinct(date, .keep_all = FALSE) %>%
          mutate(rmse_sampled = sqrt(sum((idw_cv$var1.pred - idw_cv$observed)^2)/nrow(idw_cv))) %>% as_tibble()

        #This is the RMSE value for the IDW interpolation with original points testing
        original=st_join(lcz_stations_test, st_as_sf(idw_map))
        rmse_original= data_model %>% distinct(date, .keep_all = FALSE) %>%
          mutate(rmse_original = sqrt(sum((original$IDW_airT - original$airT)^2)/nrow(original))) %>% as_tibble()
        rmse <- inner_join(rmse_original, rmse_sampled, by="date")
        return(rmse)
      }

      MapHour <- apply(ihour, 1, model_hour)
      return(MapHour)

    }

    MapDay <- pbapply::pbapply(iday, 1, model_day)
    #UHIday <- unlist(MapDay)
    return(MapDay)

  }

}

job_airT <- apply(idates, 1, AirLCZidw) #Apply the function
job_airT_list <- unlist(job_airT) #Get the raster list
job_airT_stack <- raster::stack(job_airT_list) #Or get raster stack
qtm(job_airT_list[[12]]) #plot the map









