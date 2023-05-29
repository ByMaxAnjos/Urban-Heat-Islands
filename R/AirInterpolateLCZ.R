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


#get shp ROI
shp_verify=osmdata::getbb(city, format_out = "sf_polygon", limit = 1)$polygon
if(is.null(shp_verify)) {
  study_area <- osmdata::getbb(city, format_out = "sf_polygon", limit = 1)$multipolygon# Try this first option and plot to see the city 
  study_area <- st_make_valid(study_area) %>%  st_transform(crs = mycrs)
} else {
  study_area <- osmdata::getbb(city, format_out = "sf_polygon", limit = 1)$polygon# Try this first option and plot to see the city 
  study_area <- st_make_valid(study_area) %>%  st_transform(crs = mycrs)
}
qtm(study_area) #plot map


#================================================================
#Load LCZ map 
#================================================================
# Download the LCZ global map from https://zenodo.org/record/6364594/files/lcz_filter_v1.tif?download=1
lcz_map <- raster::raster("/vsicurl/https://zenodo.org/record/6364594/files/lcz_filter_v1.tif?download=1") 

lcz_map <- raster("/Users/co2map/Documents/CO2CityMap/Berlin/Components/Building/LCZ/lcz_berlin2.tif")
lcz_map <- raster::projectRaster(lcz_map, crs = mycrs)

#Building
build <- raster::raster("/Users/co2map/Documents/CO2CityMap/Berlin/Components/building/inputs/raster/build100.tif")
names(build) <- paste0("volume")


#define the period
#month <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
imonth <- c("feb")
iyear <- c(2019)
check_day <- 2
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

AirInterpolateLCZ <- function(idates,
                           air_df = air_UCON,
                           roi = study_area,
                           spRes = 500,
                           Kriging = TRUE,
                           IDW = FALSE,
                           lcz = lcz_map,
                           isave = FALSE) {
  imonth <- idates[1]
  iyear <- idates[2]

  #Pre=processing data
  if(is.null(check_day)) {
    air_model <- air_df %>%
      openair::selectByDate(year = iyear, month = imonth) %>%
      na.omit()
  } else {
    air_model <- air_df %>%
      openair::selectByDate(year = iyear, month = imonth, day = check_day) %>%
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
        openair::selectByDate(hour = myhour) %>% as_tibble()

        #Get station points
        shp_stations <- data_model %>%
          distinct(latitude, longitude, .keep_all = T) %>%
          dplyr::select(station, latitude, longitude, airT) %>%
          sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
          st_transform(crs = mycrs) 

        #Intersect shp stations with lcz shp
        lcz_stations <- sf::st_intersection(shp_stations, lcz_shp)
       
        #Merge table to create a sf to model
        lcz_poi_mod <- inner_join(lcz_poi_get,lcz_stations %>% as_tibble() %>%
                                  dplyr::select(lcz, station, airT), by=c("lcz")) %>%
          group_by(station) %>%
          mutate(FID_station=cur_group_id()) %>%
          dplyr::select(-station, -lcz)

        # Convert LCZ map to starts
        lcz_stars <- st_as_stars(iLCZ, dimensions = "XY")
        train_mod = st_join(lcz_poi_mod, st_as_sf(lcz_stars)) %>%
          mutate(lcz = as.integer(lcz)) %>% 
          dplyr::select(airT, lcz, geometry)
        train_mod = train_mod[!is.na(train_mod$lcz),]
        st_crs(lcz_stars) <- st_crs(train_mod)

        if(Kriging==TRUE) {
          
          ## [using ordinary kriging]
          krige_vgm <- autofitVariogram(airT ~ lcz, as(train_mod, "Spatial"))$var_model
          krige_mod = gstat(formula = airT ~ lcz, model = krige_vgm$var_model, data = train_mod)
          krige_map = predict(krige_mod, newdata=lcz_stars, debug.level = 1)
          krige_map = krige_map["var1.pred",,]
          krige_map = raster::raster(rast(krige_map))
          krige_map = raster::crop(krige_map, raster::extent(roi))
          krige_map = raster::mask(krige_map, roi)
          mydate <- data_model %>% distinct(date, .keep_all = FALSE)
          mydate <- gsub("[: -]", "" , mydate$date, perl=TRUE)
          names(krige_map) <- paste0("Krige_LCZ_", mydate)
          return(krige_map)
          
        }
        
        if(IDW==TRUE) {
        
          #IDW
          idw_mod = gstat(formula = airT ~ 1, data = train_mod)
          idw_map = predict(idw_mod, lcz_stars, debug.level = 0)
          idw_map = idw_map["var1.pred",,]
          idw_map <- raster(rast(idw_map))
          mydate <- data_model %>% distinct(date, .keep_all = FALSE) %>% as.data.frame()
          mydate <- gsub("[: -]", "" , mydate$date, perl=TRUE)
          names(idw_map) <- paste0("IDW_LCZ_",mydate)
          return(idw_map)

        }
        
          # # Resample air map
          # air_resample = raster::resample(idw_map, build)
          # #raster::writeRaster(air_resample, paste0("Building/outputs/maps/UHI/", mydate, "_UHI.TIF"), format="GTiff", overwrite = TRUE)
          # #Merge building with air raster
          # air_build <- as_tibble(rasterToPoints(air_resample)) %>% set_names(c("x", "y", "airT"))
          # build_df <- as_tibble(rasterToPoints(build))
          # build_model <- inner_join(build_df, air_build, by= c("x", "y")) %>%
          #   mutate(hour = paste0(myhour))
          # 
          # #Calculate CO2 emissions
          # ECO2B <- ECO2build(build_model)
          # ECO2T_ras <- ECO2B %>% dplyr::select(x, y, ECO2_micro)
          # ECO2T_ras = raster::rasterFromXYZ(xyz = ECO2T_ras,crs = mycrs)
          # names(ECO2T_ras) <- paste0("CO2B_", mydate)
          # return(ECO2T_ras)
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

    MapHour <- pbapply::pbapply(ihour, 1, model_hour)
    UHIhour <- unlist(MapHour)
    return(UHIhour)

    }

    MapDay <- apply(iday, 1, model_day)
    UHIday <- unlist(MapDay)
    return(UHIday)

}

job_airT <- apply(idates, 1, AirInterpolateLCZ) #Apply the function
job_airT_stack <- raster::stack(unlist(job_airT)) #Or get raster stack
qtm(job_airT_stack[[c(19)]]) #plot the map









