#=======================================================================================
#
# Title:       High Spatio-temporal Resolution of Urban Heat Islands
# Author:      Dr. Max Anjos (maxanjos@campus.ul.pt)
# Description: Interpolation of air temperarure. More details on the approach are available at:
#             https://github.com/ByMaxAnjos/Urban-Heat-Islands
# Data: 2.05.2023


#' @param air_temperature_data (required). Air temperature timseries (.csv format) from different local stations that should include latitude, longitude columun
#' @param cws_data (conditionally required).
#' @return Air temperature summary, csv. table
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
  dplyr::select(date, longitude, latitude, airT) %>% 
  na.omit()

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
shp_verify=osmdata::getbb(city, format_out = "sf_polygon", limit = 1)

if(!is.null(shp_verify$geometry)) {
  study_area <- shp_verify$geometry
  study_area <- st_make_valid(study_area) %>% 
    st_transform(crs = mycrs)
} else {
  study_area <- shp_verify$multipolygon
  study_area <- st_make_valid(study_area) %>% 
    st_transform(crs = mycrs)
}

qtm(study_area) #plot map

#================================================================
#Deploy the function for interpolation of air temperature
#================================================================
'
idates is the dataframe with the defined period.
air_df is the pre-processed air temperature data
spRes is the horizontal resolution in meters
varmodel is the kriging variogram model type. It can be "Exp", "Sph", "Gau" or "Mat".
Krige if TRUE the kriging is selected.
IDW if TRUE the Inverse Distance Weighted method is selected.
Outcome this function is a data frame with the follwing metrics:n,FAC2, MB, MGE, NMB, NMGE, RMSE, r, COE, IOA.
'

#define idates
#month <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
imonth <- c("feb")
iyear <- c(2019)
check_day <- 4
idates <- expand.grid(imonth, iyear)

#Kriging

UHInterpolatEval.krige <- function(idates, air_df = air_UCON, roi = study_area, spRes = 500, 
                                   varmodel = "Sph", cv = TRUE) {
  
  imonth <- idates[1]
  iyear <- idates[2]

  air_model <- air_df %>%
    openair::selectByDate(year = iyear, month = imonth) %>% 
    na.omit()
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

  my_grid <- raster::rasterToPoints(grid_stations, spatial = TRUE)
  sp::gridded(my_grid) <- TRUE
  my_grid <- methods::as(my_grid, "SpatialPixels")
  
  #Downscale to hour
  iday <- air_model %>%
    mutate(iday = lubridate::day(date)) %>%
    distinct(iday, .keep_all = FALSE) %>%
    expand.grid()
  
  if(cv==TRUE) {
    
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
        model_cv <- distinct(data_model, longitude, latitude, .keep_all = TRUE) %>% 
          st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
          st_transform(crs = mycrs) %>% 
          dplyr::select(airT, geometry)
        krige_cv_vgm <- autofitVariogram(airT ~ 1, model_cv)
        krige_cv_mod = gstat(formula = airT ~ 1, model = krige_cv_vgm$var_model, data = model_cv)
        krige_cv = gstat.cv(krige_cv_mod, nfold = 5, verbose = FALSE)
        krige_cv = st_as_sf(krige_cv)
        
        cv_result= data_model %>% distinct(date, .keep_all = FALSE) %>%
          mutate(method = "Ordinary kriging",
                 eval = "cross validation") %>% 
          bind_cols(openair::modStats(krige_cv,  mod = "var1.pred", obs = "observed")) %>% as_tibble() %>% 
          dplyr::select(-default)
        return(cv_result)
        
      }
      MapHour <- pbapply::pbapply(ihour, 1, model_hour)
      UHIhour <- do.call(rbind.data.frame, MapHour)
      return(UHIhour)
      
    }
    
    MapDay <- apply(iday, 1, model_day)
    UHIday <- do.call(rbind.data.frame, MapDay)
    return(UHIday)
    
  }
  
  else{
    
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
        
        model_eval <- distinct(data_model, longitude, latitude, .keep_all = TRUE)
        
        #Create train and test sets
        model_split <- initial_split(model_eval)
        model_train <- training(model_split)
        model_test <- testing(model_split) 
        model_test <- st_as_sf(model_test, coords = c("longitude", "latitude"), crs = 4326) %>% 
          st_transform(crs = mycrs) %>% 
          rename(airT_actual =  airT)
        
        #Get train set
        airT_train <- model_train$airT
        get_train <- matrix(cbind(model_train$longitude,
                                  model_train$latitude),
                            length(model_train$longitude))
        
        train_df <- data.frame(x = get_train [,1], y = get_train[,2])
        sp::coordinates(train_df) = c("x", "y")
        sp::proj4string(train_df) <- sp::CRS("+proj=longlat +datum=WGS84")
        train_df <- sp::spTransform(train_df, mycrs)
        
        ## [using ordinary kriging]
        krige_mod_eval <- autofitVariogram(airT_train ~ 1, train_df)
        krige_map_eval = gstat(formula = airT_train ~ 1, model = krige_mod_eval$var_model, data = train_df)
        krige_map_eval = predict(krige_map_eval, st_as_stars(raster(my_grid)))
        krige_map_eval = krige_map_eval["var1.pred",,]
        krige_map_eval = raster::raster(rast(krige_map_eval))
        krige_map_eval = raster::crop(krige_map_eval, raster::extent(roi))
        krige_map_eval = raster::mask(krige_map_eval, roi)
        names(krige_map_eval) <- paste0("airT_pred")
        
        #Metrics
        #This is the RMSE value for the IDW interpolation with original points testing
        eval_df=st_join(model_test, st_as_sf(st_as_stars(krige_map_eval))) %>% 
          na.omit()
        eval_result= data_model %>% distinct(date, .keep_all = FALSE) %>%
          mutate(method = "Ordinary kriging") %>% 
          bind_cols(openair::modStats(eval_df,  mod = "airT_pred", obs = "airT_actual")) %>% as_tibble() %>% 
          dplyr::select(-default)
        
        return(eval_result)
        
      }
      MapHour <- pbapply::pbapply(ihour, 1, model_hour)
      UHIhour <- do.call(rbind.data.frame, MapHour)
      return(UHIhour)
      
    }
    
    MapDay <- apply(iday, 1, model_day)
    UHIday <- do.call(rbind.data.frame, MapDay)
    return(UHIday)
  }
   
}

job_krige <- pbapply::pbapply(idates, 1, UHInterpolatEval.krige)[[1]] #Apply the function 
timePlot(job_krige, pollutant = c("RMSE", "IOA", "r", "COE", "MB", "MGE", "NMGE"), y.relation = "free")



#IDW
UHInterpolatEval.IDW <- function(idates, air_df = air_UCON, roi = study_area, spRes = 500, cv=FALSE) {
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
  
  my_grid <- raster::rasterToPoints(grid_stations, spatial = TRUE)
  sp::gridded(my_grid) <- TRUE
  my_grid <- methods::as(my_grid, "SpatialPixels")
  
  #Downscale to hour
  iday <- air_model %>%
    mutate(iday = lubridate::day(date)) %>%
    distinct(iday, .keep_all = FALSE) %>%
    expand.grid()
  
  if(cv==TRUE) {
    
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
        
        model_cv <- distinct(data_model, longitude, latitude, .keep_all = TRUE) %>% 
          st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
          st_transform(crs = mycrs) %>% 
          dplyr::select(airT, geometry)
        
        #IDW
        idw_cv_mod = gstat(formula = airT ~ 1, data = model_cv)
        idw_cv = gstat.cv(idw_cv_mod, nfold = 5, verbose = FALSE)
        idw_cv = st_as_sf(idw_cv)
        cv_result= data_model %>% distinct(date, .keep_all = FALSE) %>%
          mutate(method = "IDW",
                 eval = "cross validation") %>% 
          bind_cols(openair::modStats(idw_cv,  mod = "var1.pred", obs = "observed")) %>% as_tibble() %>% 
          dplyr::select(-default)
        return(cv_result)
      }
      
      MapHour <- pbapply::pbapply(ihour, 1, model_hour)
      UHIhour <- do.call(rbind.data.frame, MapHour)
      return(UHIhour)
      
    }
    
    MapDay <- apply(iday, 1, model_day)
    UHIday <- do.call(rbind.data.frame, MapDay)
    return(UHIday)
    
  }
   
  else{
    
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
        
        model_eval <- distinct(data_model, longitude, latitude, .keep_all = TRUE)
        #Create train and test sets
        model_split <- initial_split(model_eval)
        model_train <- training(model_split) 
        model_test <- testing(model_split) 
        model_test <- st_as_sf(model_test, coords = c("longitude", "latitude"), crs = 4326) %>% 
          st_transform(crs = mycrs) %>% 
          rename(airT_actual =  airT)
        
        #Get train set
        airT_train <- model_train$airT
        get_train <- matrix(cbind(model_train$longitude,
                                  model_train$latitude),
                            length(model_train$longitude))
        
        train_df <- data.frame(x = get_train [,1], y = get_train[,2])
        sp::coordinates(train_df) = c("x", "y")
        sp::proj4string(train_df) <- sp::CRS("+proj=longlat +datum=WGS84")
        train_df <- sp::spTransform(train_df, mycrs)
        #IDW
        idw_mod = gstat(formula = airT_train ~ 1, data = train_df)
        idw_map = predict(idw_mod, st_as_stars(raster(my_grid)))
        idw_map = idw_map["var1.pred",,]
        idw_map = raster::raster(rast(idw_map))
        idw_map = raster::crop(idw_map, raster::extent(roi))
        idw_map = raster::mask(idw_map, roi)
        names(idw_map) <- paste0("airT_pred")
        
        #Metrics
        #This is the RMSE value for the IDW interpolation with original points testing
        eval_df=st_join(model_test, st_as_sf(st_as_stars(idw_map))) %>% 
          na.omit()
        eval_result= data_model %>% distinct(date, .keep_all = FALSE) %>%
          mutate(method = "IDW") %>% 
          bind_cols(openair::modStats(eval_df,  mod = "airT_pred", obs = "airT_actual")) %>% as_tibble() %>% 
          dplyr::select(-default)
        
        return(eval_result)
        
      }
      MapHour <- pbapply::pbapply(ihour, 1, model_hour)
      UHIhour <- do.call(rbind.data.frame, MapHour)
      return(UHIhour)
      
    }
    
    MapDay <- apply(iday, 1, model_day)
    UHIday <- do.call(rbind.data.frame, MapDay)
    return(UHIday)
  }
  
}

job_idw <- apply(idates, 1, UHInterpolatEval.IDW)[[1]]

timePlot(job_idw, pollutant = c("RMSE", "IOA", "r", "COE", "MB", "MGE", "NMGE"), y.relation = "free")


# Plot
tm_shape(krige_map) + 
  tm_raster(n=10, palette = "RdBu", auto.palette.mapping = FALSE,
            title="Predicted Air temperature \n(in degree)") + 
  tm_shape(model_test) + tm_dots(size=0.2) +
  tm_text("airT_actual", just="left", xmod=.5, size = 0.7) +
  tm_legend(legend.outside=TRUE)



