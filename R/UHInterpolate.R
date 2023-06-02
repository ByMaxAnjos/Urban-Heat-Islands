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
library(tidymodels)
library(stars)
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

#define the period
#month <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
imonth <- c("jul")
iyear <- c(2018)
check_day <- 3
idates <- expand.grid(imonth, iyear)


UHInterpolate <- function(idates,
                           air_df = air_UCON,
                           roi = study_area,
                           spRes = 500,
                           kriging = FALSE,
                           varmodel = "Sph",
                           IDW = TRUE,
                           isave = TRUE) {
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
     
      if(kriging==TRUE) {
        
        ## [using ordinary kriging]
        krige_mod <- autofitVariogram(airT ~ 1, get_coord)
        krige_mod = gstat(formula = airT ~ 1, model = krige_mod$var_model, data =get_coord)
        krige_map = predict(krige_mod, st_as_stars(raster(my_grid)))
        krige_map = krige_map["var1.pred",,]
        krige_map = raster::raster(rast(krige_map))
        krige_map = raster::crop(krige_map, raster::extent(roi))
        krige_map = raster::mask(krige_map, roi)
        mydate <- data_model %>% distinct(date, .keep_all = FALSE) %>% as.data.frame()
        mydate <- gsub("[: -]", "" , mydate$date, perl=TRUE)
        names(krige_map) <- paste0("Krige_", mydate)
        return(krige_map)
      
      }
      
      if(isave==TRUE & kriging==TRUE) {
        
        # Create a folder name using paste0
        folder <- paste0("UHI_Output/")
        
        # Check if the folder exists
        if (!dir.exists(folder)) {
          # Create the folder if it does not exist
          dir.create(folder)
        }
        
        file <- paste0(folder,iyear,imonth,myday,myhour,"Krige.TIF")
        raster::writeRaster(krige_map,file, format="GTiff", overwrite = TRUE)
        
      }
      
      
      if(IDW==TRUE){
        
        #IDW
        idw_mod = gstat(formula = airT ~ 1, data = get_coord)
        idw_map = predict(idw_mod, st_as_stars(raster(my_grid)))
        idw_map = idw_map["var1.pred",,]
        idw_map = raster::raster(rast(idw_map))
        idw_map = raster::crop(idw_map, raster::extent(roi))
        idw_map = raster::mask(idw_map, roi)
        mydate <- data_model %>% distinct(date, .keep_all = FALSE) %>% as.data.frame()
        mydate <- gsub("[: -]", "" , mydate$date, perl=TRUE)
        names(idw_map) <- paste0("IDW_", mydate)
        return(idw_map)
        
      }
      
      if(isave == TRUE & IDW == TRUE) {
        
        # Create a folder name using paste0
        folder <- paste0("UHI_Output/")
        
        # Check if the folder exists
        if (!dir.exists(folder)) {
          # Create the folder if it does not exist
          dir.create(folder)
        }
        
        file <- paste0(folder,iyear,imonth,myday,myhour,"IDW.TIF")
        raster::writeRaster(idw_map,file, format="GTiff", overwrite = TRUE)
        
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

job_airT <- apply(idates, 1, UHInterpolate) #Apply the function 
job_airT_stack <- raster::stack(unlist(job_airT)) #Or get raster stack
qtm(job_airT_stack[[7]]) #plot the map

