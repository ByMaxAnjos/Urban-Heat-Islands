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
#Define your path
#setwd("Path/")

#================================================================
#Load air temperature data and pre-processing
#================================================================

#CWS data
air_cws_2018 <- fread("/Users/co2map/Documents/CO2CityMap/Berlin/Components/Building/inputs/data/cws_berlin_ta_level_o1_2018_July.csv")
#air_cws_2019 <- fread("/Users/co2map/Documents/CO2CityMap/Berlin/Components/Building/inputs/data/cws_berlin_ta_level_o1_2019_July.csv")

#Pre=processing data
air_cws_2018 <- air_cws_2018 %>%
  rename(date=time,
         longitude = lon,
         latitude = lat,
         airT = ta_int) %>%
  dplyr::select(date, longitude, latitude, airT) 

#UCON data
air_UCON <- fread("/Users/co2map/Documents/CO2CityMap/Berlin/Components/building/inputs/data/airT_UCON_2015_2022.csv") %>%
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


#define the period
#month <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
imonth <- c("jul")
iyear <- c(2018)
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

AirInterpolate <- function(idates,
                           air_df = air_cws_2018,
                           roi = study_area,
                           spRes = 500,
                           impute = FALSE,
                           varmodel = "Sph",
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
  air_model <- air_df %>%
    openair::selectByDate(year = iyear, month = imonth) %>% 
    na.omit()
  
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
        openair::selectByDate(hour = myhour)
      
      airT <- data_model$airT
      get_coord <- matrix(cbind(data_model$longitude,
                                data_model$latitude),
                          length(data_model$longitude))
      
      get_coord_df <- data.frame(x = get_coord [,1], y = get_coord[,2])
      sp::coordinates(get_coord_df) = c("x", "y")
      sp::proj4string(get_coord_df) <- sp::CRS("+proj=longlat +datum=WGS84")
      get_coord <- sp::spTransform(get_coord_df, mycrs)
      
      grid_stations <- raster::raster(raster::extent(roi), res = spRes)
      raster::values(grid_stations) <- 1:raster::ncell(grid_stations)
      raster::crs(grid_stations) <- mycrs
      
      # Convert to spatial pixel
      my_grid <- raster::rasterToPoints(grid_stations, spatial = TRUE)
      sp::gridded(my_grid) <- TRUE
      my_grid <- methods::as(my_grid, "SpatialPixels")
      
      vgm_air = gstat::variogram(object = airT ~ 1, data = get_coord)# set a variogram
      fit_var_air = gstat::fit.variogram(object = vgm_air, gstat::vgm(varmodel)) # fit a variogram
      fit_var_air$range[fit_var_air$range < 0] <- abs(fit_var_air$range)[2]
      kriging_air = gstat::krige(airT ~ 1, locations = get_coord, newdata = my_grid,
                                 model = fit_var_air, debug.level = 0)
      kriging_air = raster::raster(kriging_air)
      kriging_air = raster::crop(kriging_air, raster::extent(roi))
      kriging_air = raster::mask(kriging_air, roi)
      mydate <- data_model %>% distinct(date, .keep_all = FALSE) %>% as.data.frame()
      mydate <- gsub("[: -]", "" , mydate$date, perl=TRUE)
      names(kriging_air) <- paste0("airT_",mydate)

      return(kriging_air)
      
      if(isave == TRUE) {
        
        # Create a folder name using paste0
        folder <- paste0("ZCCM_Output/")
        
        # Check if the folder exists
        if (!dir.exists(folder)) {
          # Create the folder if it does not exist
          dir.create(folder)
        }
        
        file <- paste0(folder,iyear,imonth,myday,myhour,"AirT.TIF")
        raster::writeRaster(kriging_air,file, format="GTiff", overwrite = TRUE)
        
      }
      
    }
    
    MapHour <- pbapply::pbapply(ihour, 1, model_hour)
    return(MapHour)
    
  }
  
  MapDay <- apply(iday, 1, model_day)
  UHIday <- unlist(MapDay)
  return(MapDay)
  
}

job_airT <- apply(idates, 1, AirInterpolate) #Apply the function 
job_airT_list <- unlist(job_airT) #Get the raster list
job_airT_stack <- raster::stack(job_airT_list) #Or get raster stack
qtm(job_airT_list[[7]]) #plot the map



