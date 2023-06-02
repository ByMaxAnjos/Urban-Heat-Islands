#=======================================================================================
#
# Title:       High Spatio-temporal Resolution of Urban Heat Islands
# Author:      Dr. Max Anjos (maxanjos@campus.ul.pt)
# Description: Interpolation of air temperarure. More details on the approach are available at:
#             https://github.com/ByMaxAnjos/Urban-Heat-Islands
# Data: 27.05.2023


#' @param air_temperature_data (required). Air temperature timseries (.csv format) from different local stations that should include latitude, longitude column
#' @param LCZ_map (required).
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
library(terra)
library(tidymodels)
library(stars)
#Define your path
#setwd("Path/")


#================================================================
#FUNCTION X: get your LCZ map 
#================================================================
'
Function to download and process Local Climate Zone Classification data from a specified region of interest (ROI)
Input: city (optional) - name of the city to download data for
roi (optional) - polygon object defining the region of interest
Output: raster object with the LCZ classes, with the name "LCZ_<city/ROI>"
'

getLCZmap <- function(city=NULL, roi = NULL) {
    # Validate inputs
    if (is.null(city) & is.null(roi)) {
      stop("Error: provide either a city name or a roi polygon")
    } else if (!is.null(city) & !is.character(city)) {
      stop("Error: city input must be a character string")
    } else if (!is.null(roi) & !inherits(roi, "sf")) {
      stop("Error: ROI input must be a polygon object of class sf")
    }
    
    if(!is.null(city)) {
      
      # Get study area polygon from OpenStreetMap data
      shp_verify <- osmdata::getbb(city, format_out = "sf_polygon", limit = 1, featuretype = "city")
      
      # Check if polygon was obtained successfully
      if(!is.null(shp_verify$geometry) & !inherits(shp_verify, "list")) {
        study_area <- shp_verify$geometry
        study_area <- st_make_valid(study_area) %>%
          st_as_sf() %>% 
          st_transform(crs = "+proj=longlat +datum=WGS84 +no_defs")
      } else {
        study_area <- shp_verify$multipolygon
        study_area <- st_make_valid(study_area) %>%
          st_as_sf() %>%
          st_transform(crs="+proj=longlat +datum=WGS84 +no_defs")
      }
      # Download the LCZ global map from https://zenodo.org/record/6364594/files/lcz_filter_v1.tif?download=1
    #lcz_download <- raster::raster("/vsicurl/https://zenodo.org/record/6364594/files/lcz_filter_v1.tif?download=1")
    lcz_ras <- raster::crop(lcz_download, raster::extent(study_area))
    lcz_ras <- raster::mask(lcz_ras, study_area)
    names(lcz_ras) <- paste0("LCZ_",city)
    return(lcz_ras)
  } else {
    # Download the LCZ global map from https://zenodo.org/record/6364594/files/lcz_filter_v1.tif?download=1
    #lcz_download <- raster::raster("/vsicurl/https://zenodo.org/record/6364594/files/lcz_filter_v1.tif?download=1") 
    roi_crs <- roi %>% st_as_sf() %>% st_transform(crs = "+proj=longlat +datum=WGS84 +no_defs")
    lcz_ras <- raster::crop(lcz_download, raster::extent(roi_crs))
    lcz_ras <- raster::mask(lcz_ras, roi_crs)
    names(lcz_ras) <- "LCZ_ROI"
    return(lcz_ras)
    
  }
  
}

lcz_ber <- getLCZmap(city="Berlin")

#================================================================
#FUNCTION X: plot your LCZ map 
#================================================================


plotLCZ <- function(x) {
  
  lcz_map <- x
  
  lczClass <- ratify(lcz_map)
  
  rat <- levels(lczClass)[[1]]
  
  ID <- c(seq(1, 10, 1), seq(11, 17)) %>% as_tibble() %>% set_names("ID")
  
  lcz.name <- c("Compact highrise", "Compact midrise", "Compact lowrise", "Open highrise",
                "Open midrise", "Open lowrise", "Lightweight low-rise", "Large lowrise",
                "Sparsely built", "Heavy Industry", "Dense trees", "Scattered trees",
                "Bush, scrub", "Low plants", "Bare rock or paved", "Bare soil or sand", "Water") %>% as_tibble() %>% set_names("lcz.name")
  lcz.col <- c("#910613", "#D9081C", "#FF0A22", "#C54F1E", "#FF6628", "#FF985E",
               "#FDED3F", "#BBBBBB", "#FFCBAB", "#565656", "#006A18", "#00A926",
               "#628432", "#B5DA7F", "#000000", "#FCF7B1", "#656BFA") %>% as_tibble() %>% set_names("lcz.col")
  
  lcz_df <- bind_cols(ID, lcz.name, lcz.col) %>% 
    inner_join(rat, by = "ID")
  
  rat$classes <- c(lcz_df$lcz.name)
  levels(lczClass) <- rat
  qualPal <- c(lcz_df$lcz.col)
  qualTheme <- rasterVis::rasterTheme(region = qualPal,
                                      panel.background = list(col = 'white')
  )

  # rasterVis::levelplot(lczClass, at = seq(min(lcz_df$ID),max(lcz_df$ID)), 
  #                      par.settings = qualTheme)
  
  tm_shape(lczClass) +
    tm_raster(palette = qualPal, labels = lcz_df$lcz.name, title = "LCZ") +
    tm_layout(legend.show = TRUE,
              title.color = "#3f1651", title.size = 1.5,
              frame = TRUE, legend.outside = TRUE, legend.text.size = 1) +
    #tm_scale_bar(position = "right", lwd = 1, color.dark = "black", color.light = "white", text.size = 0.5)+
    tm_credits("Source: ©ZoomCityCarbonModel,https://github.com/ByMaxAnjos/Urban-Heat-Islands\nData:Demuzere et al.(2022)/https://doi.org/10.5194/essd-14-3835-2022\nOpenStreetMap® contributions, 2023",
               size = 0.4, position=c("left","bottom"), col = "#3f1651", bg.alpha = 0.5)+
    #tm_compass(type="arrow", position=c("right", "top"), show.labels = 1) +
    tm_graticules(lwd = .1) +
    qtm(study_area, fill = NULL)
  #tmap_save(iplot, filename = "fig-NEEmap.png", asp = 0)
  

}

g <- plotLCZ(lcz_ber)

#Salve map
tmap_save(g, "lcz_city.png", asp = 0)


#================================================================
#FUNCTION X: get parameters from LCZ
#================================================================
'
This nice function gets all LCZ parameters from Steawrt and Oke (2012)
and convert them to shpafile or, if ras=TRUE,to raster stack.
'
getLCZparameters <- function(lcz_map, ras = FALSE) {
  #lcz.id <- c(seq(1, 10, 1), seq(101, 107))
  lcz <- c(seq(1, 10, 1), seq(11, 17))
  lcz.code <- c(seq(1, 10, 1), "A", "B", "C", "D", "E", "F", "G")
  # lcz.name <- c('compact_high-rise', 'compact_midrise', 'compact_low-rise',
  #               'open_high-rise', 'open_midrise', 'open_low-rise',
  #               'lightweight_low-rise', 'large_low-rise', 'sparsely_built',
  #               'heavy_industry', 'dense_trees', 'scattered_trees', 'bush_scrub',
  #               'low_plants', 'bare_rock_paved', 'bare_soil_sand', 'water')
  
  lcz.name <- c("Compact highrise", "Compact midrise", "Compact lowrise", "Open highrise",
                "Open midrise", "Open lowrise", "Lightweight low-rise", "Large lowrise",
                "Sparsely built", "Heavy Industry", "Dense trees", "Scattered trees",
                "Bush, scrub", "Low plants", "Bare rock or paved", "Bare soil or sand", "Water")
  lcz.col <- c("#910613", "#D9081C", "#FF0A22", "#C54F1E", "#FF6628", "#FF985E",
               "#FDED3F", "#BBBBBB", "#FFCBAB", "#565656", "#006A18", "#00A926",
               "#628432", "#B5DA7F", "#000000", "#FCF7B1", "#656BFA")
  
  #LCZ parameters
  SVF.min <- c(0.2, 0.3, 0.2, 0.5, 0.5, 0.6, 0.2, 0.75, 0.85, 0.6, 0.35, 0.5, 0.7, rep(0.9, 4))
  SVF.max <- c(0.4, 0.6, 0.6, 0.7, 0.8, 0.9, 0.5, 0.75, 0.85, 0.9, 0.35, 0.8, 0.9, rep(0.9, 4))
  aspect.ratio.min <- c(3, 0.75, 0.75, 0.75, 0.3, 0.3, 1, 0.1, 0.1, 0.2, 1.5, 0.25, 0.25, rep(0.1, 4))
  aspect.ratio.max <- c(3, 2, 1.5, 1.25, 0.75, 0.75, 2, 0.3, 0.25, 0.5, 1.5, 0.75, 1.0, rep(0.1, 4))
  build.frac.min <- c(40, 40, 40, rep(20,3), 60, 30, 10, 20, rep(9, 7))
  build.frac.max <- c(60, 70, 70, rep(40,3), 90, 50, 20, 30, rep(9, 7))
  imp.frac.min <- c(40, 40, 40, rep(20, 3), 60, 30, 10, 20, rep(0, 7))
  imp.frac.max <- c(60, 70, 70, rep(40, 3), 90, 50, 20, 30, rep(10, 7))
  veg.frac.max <- c(10, 20, 30, 40, 40, 60, 30, 20, 80, 50, rep(100, 4), 10, 100, 100)
  veg.frac.min <- c(0, 0, 0, 30, 20, 30, 0, 0, 60, 40, 90, 90, 90, 90, 0, 90, 90)
  tree.frac.min <- c(rep(0, 10), 90, 90, rep(0, 5))
  tree.frac.max <- c(rep(0, 10), 100, 100, rep(0, 5))
  height.roug.min <- c(26, 10, 3, 26, 10, 3, 2, 3, 3, 5, 3, 3, 2.9, 0.9, 0.24, 0.23,  0)
  height.roug.max <- c(26, 25, 10, 26, 25, 10, 4, 10, 10, 15, 30, 15, 2.9, 0.9, 0.24, 0.23, 0)
  terra.roug.min <- c(8, 6, 6, 7, 5, 5, 4, 5, 5, 5, 8, 5, 4, 3, 1, 1, 1)
  terra.roug.max <- c(8, 7, 6, 8, 6, 6, 5, 5, 6, 6, 8, 6, 5, 4, 2, 2, 1)
  surf.admit.min <- c(1.500, 1.500, 1.200, 1.400, 1.400, 1.200, 800, 1.200, 1.000, 1.000, 0, 1.000, 700, 1.200, 1.200, 600, 1.500)
  surf.admit.max <- c(1.800, 2.000, 1.800, 1.800, 2.000, 1.800, 1.500, 1.800, 1.800, 2.5000, 0, 1.800, 1.500, 1.600, 2.500, 1.400, 1.500)
  surf.albedo.min <- c(rep(0.10, 3), rep(0.12, 3), rep(0.15, 2), rep(0.12, 2), 0.10, rep(0.15, 4), 0.20, 0.02)
  surf.albedo.max <- c(rep(0.20, 3), rep(0.25, 3), 0.35, 0.25, 0.25, 0.20, 0.20, 0.25, 0.30, 0.25, 0.30, 0.35, 0.10)
  antrop.heat.min <- c(50, 74, 74, 49, 24, 24, 34, 49, 9, 310, rep(0, 7))
  antrop.heat.max <- c(300, 74, 74, 49, 24, 24, 34, 49, 9, 310, rep(0, 7))
  
  # lcz.col <- c('#8c0000', '#d10000', '#ff0100', '#be4d01', '#ff6602', '#ff9955',
  #              '#faee05', '#bcbcbc', '#ffccaa', '#555555', '#006a01', '#01aa00',
  #              '#648526', '#b9db79', '#000000', '#fbf7ae', '#6a6aff')
  lcz.df <- data.frame(lcz, lcz.name, lcz.code, lcz.col, SVF.min, SVF.max, aspect.ratio.min, aspect.ratio.max, build.frac.min, build.frac.max,
                       imp.frac.min, imp.frac.max, veg.frac.min, veg.frac.max, tree.frac.min, tree.frac.max,
                       height.roug.min, height.roug.max, terra.roug.min, terra.roug.max, surf.admit.min, surf.admit.max, surf.albedo.min, surf.albedo.max,
                       antrop.heat.min, antrop.heat.max,
                       stringsAsFactors = F) %>%
    mutate(z0 = ifelse(lcz.code %in% c("G"), 0.0002, #Get z0
                       ifelse(lcz.code %in% c("E", "F"), 0.0005,
                              ifelse(lcz.code=="D", 0.03,
                                     ifelse(lcz.code %in% c(7, "C"), 0.10,
                                            ifelse(lcz.code %in% c(8, "B"), 0.25,
                                                   ifelse(lcz.code %in% c(2, 3, 5, 6, 9, 10), 0.5,
                                                          ifelse(lcz.code %in% c(2, 4), 1.0,
                                                                 ifelse(lcz.code %in% c(1, "A"), 2, ""))))))))) %>%
    mutate(SVF.mean = round((SVF.min + SVF.max)/2, digits = 2),
           aspect.ratio.mean = (aspect.ratio.min + aspect.ratio.max)/2,
           build.frac.mea = (build.frac.min + build.frac.max)/2,
           imp.frac.mean = (imp.frac.min + imp.frac.max)/2,
           veg.frac.mean = (veg.frac.min + veg.frac.max)/2,
           tree.frac.mean = (tree.frac.min +tree.frac.max)/2,
           height.roug.mean = (height.roug.min + height.roug.max)/2,
           terra.roug.mean = (terra.roug.min + terra.roug.max)/2,
           surf.admit.mean = (surf.admit.min + surf.admit.max)/2,
           surf.albedo.mean = (surf.albedo.min + surf.albedo.max)/2,
           antrop.heat.mean = (antrop.heat.min + antrop.heat.max)/2
    )
  #Preprocessing raster
  names(lcz_map) <- "lcz"
  lcz_shp <- terra::as.polygons(rast(lcz_map)) %>% st_as_sf()
  lcz_result <- inner_join(lcz_shp, lcz.df, by="lcz")
  
  if(ras==TRUE){
    ras_map <- pbapply::pblapply(6:ncol(lcz_result)-1, FUN=function(i)
      stars::st_rasterize(lcz_result[,i]) %>% 
        terra::rast() %>% raster::raster() %>% 
        terra::resample(lcz_map)
      )
    ras_stack <- raster::stack(ras_map)
    names(ras_stack) <- colnames(lcz_result)[5:38]
    return(ras_stack)
  } else {
    return(lcz_result) 
  }
}
myLCZpar <- getLCZparameters(lcz_ara, ras = TRUE)

# Plot
tm_shape(myLCZpar) + 
  tm_raster(n=10, palette = "RdBu",
            title="LCZ parameters") + 
  tm_layout(legend.show = FALSE)

#================================================================
#FUNCTION 2: calculate the time series of anomaly between LCZ
#================================================================
'
This function calculate the anomaly between the classifed LCZ.
It uses the same dataframe of air temperature and station points. 
'




#define idates
#month <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
imonth <- c("feb")
iyear <- c(2019)
check_day <- 4
idates <- expand.grid(imonth, iyear)




#================================================================
#FUNCTION X: Interpolation of air temperature based on LCZ map
#================================================================
'
idates is the dataframe with the defined period.
air_df is the pre-processed air temperature data
spRes is the horizontal resolution in meters
impute if TRUE impute the missing air "NaN" crws data.
varmodel is the variogram model type. It can be "Exp", "Sph", "Gau" or "Mat".
isave if TRUE creates a new folder "ZCCM_Output" in you path and saves the rasters
'

#Temperature data
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

#Get lcz map
lcz_map <- raster("/Users/co2map/Documents/CO2CityMap/Berlin/Components/Building/LCZ/lcz_berlin2.tif")
lcz_map <- raster::projectRaster(lcz_map, crs = mycrs)
qtm(lcz_map)

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
        mydate <- data_model %>% distinct(date, .keep_all = FALSE) %>% as.data.frame()
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


#================================================================
#FUNCTION X: Evaluate the LCZ - Interpolate air temperature
#================================================================
'
This function evaluates the combination of LCZ-interpolate
It uses the same dataframe of air temperature and LCZ map.
'

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


#define idates
#month <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
imonth <- c("feb")
iyear <- c(2019)
check_day <- 4
idates <- expand.grid(imonth, iyear)

#LCZ
UHInterpolatEval.LCZ <- function(idates,
                                 lcz = lcz_map,
                                 air_df = air_UCON,
                                 roi = study_area,
                                 spRes = 500,
                                 Kriging.eval = FALSE,
                                 kriging.cv = FALSE,
                                 IDW.eval = FALSE,
                                 IDW.cv = TRUE) {
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
    #raster::mask(roi) %>%
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
  
  if(Kriging.eval==TRUE) {
    
    
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
        lcz_stations_train <- sf::st_intersection(shp_stations, lcz_shp)
        qtm(lcz_shp) + qtm(shp_stations)
        model_split <- initial_split(lcz_stations_train, prop = .7, breaks = 2)
        model_train <- training(model_split)
        model_test <- testing(model_split) %>% st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
          st_transform(crs = mycrs) %>% 
          rename(airT_actual =  airT)
        
        
        #Merge table to create a sf to model
        lcz_poi_mod <- inner_join(lcz_poi_get, model_train %>% as_tibble() %>%
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
        
        ## [using ordinary kriging]
        krige_vgm <- autofitVariogram(airT ~ lcz, as(train_mod, "Spatial"))$var_model
        krige_mod = gstat(formula = airT ~ lcz, model = krige_vgm$var_model, data = train_mod)
        krige_map = predict(krige_mod, newdata=lcz_stars, debug.level = 1)
        krige_map = krige_map["var1.pred",,]
        krige_map = raster::raster(rast(krige_map))
        # krige_map = raster::crop(krige_map, raster::extent(roi))
        # krige_map = raster::mask(krige_map, roi)
        mydate <- data_model %>% distinct(date, .keep_all = FALSE) %>% as.data.frame()
        mydate <- gsub("[: -]", "" , mydate$date, perl=TRUE)
        names(krige_map) <- paste0("airT_pred")
        
        
        #Metrics
        #This is the RMSE value for the IDW interpolation with original points testing
        eval_df=st_join(model_test, st_as_sf(st_as_stars(krige_map))) %>% 
          na.omit() %>% as_tibble() %>% dplyr::select(-geometry)
        eval_result= data_model %>% distinct(date, .keep_all = FALSE) %>%
          mutate(method = "LCZ-Ordinary kriging") %>% 
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
  
  if(kriging.cv==TRUE) {
    
    
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
        
        model_split <- initial_split(shp_stations)
        model_train <- training(model_split)
        model_test <- testing(model_split) %>% st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
          st_transform(crs = mycrs) %>% 
          rename(airT_actual =  airT)
        
        #Intersect shp stations with lcz shp
        lcz_stations_train <- sf::st_intersection(model_train, lcz_shp)
        
        #Merge table to create a sf to model
        lcz_poi_mod <- inner_join(lcz_poi_get,lcz_stations_train %>% as_tibble() %>%
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
        
        ## [using ordinary kriging]
        lcz_vgm <- autofitVariogram(airT ~ lcz, as(train_mod, "Spatial"))
        lcz_mod = gstat(formula = airT ~ lcz, model = lcz_vgm$var_model, data = lcz_stations_train)
        lcz_cv = gstat.cv(lcz_mod, nfold = 5, verbose = FALSE)
        lcz_cv = st_as_sf(lcz_cv)
        
        cv_result= data_model %>% distinct(date, .keep_all = FALSE) %>%
          mutate(method = "LCZ-Universal kriging",
                 eval = "cross validation") %>% 
          bind_cols(openair::modStats(lcz_cv,  mod = "var1.pred", obs = "observed")) %>% as_tibble() %>% 
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
  
  if(IDW.eval==TRUE) {
    
    
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
        
        model_split <- initial_split(shp_stations)
        model_train <- training(model_split)
        model_test <- testing(model_split) %>% st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
          st_transform(crs = mycrs) %>% 
          rename(airT_actual =  airT)
        
        #Intersect shp stations with lcz shp
        lcz_stations_train <- sf::st_intersection(model_train, lcz_shp)
        
        #Merge table to create a sf to model
        lcz_poi_mod <- inner_join(lcz_poi_get,lcz_stations_train %>% as_tibble() %>%
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
        
        #IDW
        idw_mod = gstat(formula = airT ~ 1, data = train_mod)
        idw_map = predict(idw_mod, lcz_stars, debug.level = 0)
        idw_map = idw_map["var1.pred",,]
        idw_map <- raster(rast(idw_map))
        mydate <- data_model %>% distinct(date, .keep_all = FALSE) %>% as.data.frame()
        mydate <- gsub("[: -]", "" , mydate$date, perl=TRUE)
        names(idw_map) <- paste0("airT_pred")
        
        #Metrics
        #This is the RMSE value for the IDW interpolation with original points testing
        eval_df=st_join(model_test, st_as_sf(st_as_stars(idw_map))) %>% 
          na.omit() %>% as_tibble() %>% dplyr::select(-geometry)
        eval_result= data_model %>% distinct(date, .keep_all = FALSE) %>%
          mutate(method = "LCZ-IDW") %>% 
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
  
  if(IDW.cv==TRUE) {
    
    
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
        
        model_split <- initial_split(shp_stations)
        model_train <- training(model_split)
        model_test <- testing(model_split) %>% st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
          st_transform(crs = mycrs) %>% 
          rename(airT_actual =  airT)
        
        #Intersect shp stations with lcz shp
        lcz_stations_train <- sf::st_intersection(model_train, lcz_shp)
        
        #Merge table to create a sf to model
        lcz_poi_mod <- inner_join(lcz_poi_get,lcz_stations_train %>% as_tibble() %>%
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
        
        #IDW
        idw_cv_mod = gstat(formula = airT ~ 1, data = train_mod)
        idw_cv = gstat.cv(idw_cv_mod, nfold = 5, verbose = FALSE)
        idw_cv = st_as_sf(idw_cv)
        cv_result= data_model %>% distinct(date, .keep_all = FALSE) %>%
          mutate(method = "LCZ-IDW",
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
  
}

job_lcz <- apply(idates, 1, UHInterpolatEval.LCZ)[[1]] #Apply the function

#Plot the metrics
timePlot(job_airT, pollutant = c("RMSE", "IOA", "r", "COE", "MB", "MGE", "NMGE"), y.relation = "free")


#================================================================
#FUNCTION X: Estimate the CO2 emissions from buildings using LCZ map
#================================================================
'
idates is the dataframe with the defined period.
air_df is the pre-processed air temperature data
spRes is the horizontal resolution in meters
impute if TRUE impute the missing air "NaN" crws data.
varmodel is the variogram model type. It can be "Exp", "Sph", "Gau" or "Mat".
isave if TRUE creates a new folder "ZCCM_Output" in you path and saves the rasters
'

#Building
build <- raster::raster("/Users/co2map/Documents/CO2CityMap/Berlin/Components/building/inputs/raster/build100.tif")


CO2BuildLCZ <- function(idates,
                              air_df = air_UCON,
                              lcz = lcz_map,
                              build = build,
                              roi = study_area,
                              spRes = 100,
                              Kriging = FALSE,
                              IDW = TRUE) {
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
  
  if(kriging==TRUE) {
    
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
   
          ## [using ordinary kriging]
          krige_vgm <- autofitVariogram(airT ~ lcz, as(train_mod, "Spatial"))$var_model
          krige_mod = gstat(formula = airT ~ lcz, model = krige_vgm$var_model, data = train_mod)
          krige_map = predict(krige_mod, newdata=lcz_stars, debug.level = 1)
          krige_map = krige_map["var1.pred",,]
          krige_map = raster::raster(rast(krige_map))
          krige_map = raster::crop(krige_map, raster::extent(roi))
          krige_map = raster::mask(krige_map, roi)
          mydate <- data_model %>% distinct(date, .keep_all = FALSE) %>% as.data.frame()
          mydate <- gsub("[: -]", "" , mydate$date, perl=TRUE)
          names(krige_map) <- paste0("Krige_LCZ_", mydate)
          
          # Resample air map
          air_resample = raster::resample(krige_map, build)
          #Rename raster to "volume"
          names(build) <- paste0("volume")
          #Merge building with air raster
          air_build <- as_tibble(rasterToPoints(air_resample)) %>% set_names(c("x", "y", "airT"))
          build_df <- as_tibble(rasterToPoints(build))
          build_model <- inner_join(build_df, air_build, by= c("x", "y")) %>%
            mutate(hour = paste0(myhour))
          #Formula to calculate building CO2 emissions
          ECO2build <- function(x) {
            x <- mutate(x, ECO2_micro = as.numeric(
              ifelse(hour == 0, abs(-1.57*max(0, 9.9 - airT)-3.6/volume), 
                     ifelse(hour == 1, abs(-1.73*max(0, 12.8 - airT)- 4.2/volume),
                            ifelse(hour == 2, abs(-1.73*max(0,12 - airT)-4.64/volume),
                                   ifelse(hour == 3, abs(-1.5*max(0, 13.1 - airT)-3.44/volume), 
                                          ifelse(hour == 4, abs(-1.38*max(0, 12.4 - airT)-3.86/volume), 
                                                 ifelse(hour == 5, abs(-1.7*max(0, 12.9 - airT)-2.46/volume),
                                                        ifelse(hour == 6, abs(-2.55*max(0, 12.7 - airT)-0.77/volume), 
                                                               ifelse(hour == 7, abs(-2.29*max(0, 10.2 - airT)-2.18/volume), 
                                                                      ifelse(hour == 8, abs(-2.88*max(0, 12.7 - airT)- (-1.5)/volume),
                                                                             ifelse(hour == 9, abs(-2.97*max(0, 12.3 - airT)-(-1.8)/volume), 
                                                                                    ifelse(hour == 10, abs(-1.27*max(0, 13.3 - airT)-0.59/volume), 
                                                                                           ifelse(hour == 11, abs(-1.45*max(0, 15.3 - airT)-0.57/volume),
                                                                                                  ifelse(hour == 12, abs(-1.48*max(0, 15 - airT)-1.1/volume), 
                                                                                                         ifelse(hour == 13, abs(-1.45*max(0, 13.7 - airT)-1.7/volume), 
                                                                                                                ifelse(hour == 14, abs(-1.34*max(0, 13.2 - airT)-2.11/volume),
                                                                                                                       ifelse(hour == 15, abs(-1.35*max(0, 13.5 - airT)-3.22/volume), 
                                                                                                                              ifelse(hour == 16, abs(-1.36*max(0, 13 - airT)-4.39/volume), 
                                                                                                                                     ifelse(hour == 17, abs(-1.29*max(0, 12 - airT)-5.4/volume),
                                                                                                                                            ifelse(hour == 18, abs(-1.44*max(0, 10.6 - airT)-7.3/volume), 
                                                                                                                                                   ifelse(hour == 19, abs(-1.44*max(0, 7 - airT)-9.77/volume), 
                                                                                                                                                          ifelse(hour == 20, abs(-2.7*max(0, 6 - airT)-6.87/volume),
                                                                                                                                                                 ifelse(hour == 21, abs(-1.24*max(0, 8 - airT)-4.29/volume), 
                                                                                                                                                                        ifelse(hour == 22, abs(-0.48*max(0, 10.3 - airT)-3.68/volume),
                                                                                                                                                                               ifelse(hour == 23, abs(-1.69*max(0, 12.4 - airT)-3.6/volume), ""))))))))))))))))))))))))))
            
            return(x)
          }
          
          #Calculate CO2 emissions
          ECO2B <- ECO2build(build_model)
          ECO2B_ras <- ECO2B %>% dplyr::select(x, y, ECO2_micro)
          ECO2B_ras = raster::rasterFromXYZ(xyz = ECO2B_ras,crs = mycrs)
          names(ECO2B_ras) <- paste0("CO2B_", mydate)
          qtm(ECO2B_ras)
          return(ECO2B_ras)
      }
      
      MapHour <- pbapply::pbapply(ihour, 1, model_hour)
      UHIhour <- unlist(MapHour)
      return(UHIhour)
      
    }
    
    MapDay <- apply(iday, 1, model_day)
    UHIday <- unlist(MapDay)
    return(UHIday)
  }
  
  if(IDW==TRUE) {
    
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
          
          #IDW
          idw_mod = gstat(formula = airT ~ 1, data = train_mod)
          idw_map = predict(idw_mod, lcz_stars, debug.level = 0)
          idw_map = idw_map["var1.pred",,]
          idw_map <- raster(rast(idw_map))
          mydate <- data_model %>% distinct(date, .keep_all = FALSE) %>% as.data.frame()
          mydate <- gsub("[: -]", "" , mydate$date, perl=TRUE)
          names(idw_map) <- paste0("IDW_LCZ_",mydate)
          
          # Resample air map
          air_resample = raster::resample(idw_map, build)
          #Rename raster to "volume"
          names(build) <- paste0("volume")
          #Merge building with air raster
          air_build <- as_tibble(rasterToPoints(air_resample)) %>% set_names(c("x", "y", "airT"))
          build_df <- as_tibble(rasterToPoints(build))
          build_model <- inner_join(build_df, air_build, by= c("x", "y")) %>%
            mutate(hour = paste0(myhour))
          
          #Formula to calculate building CO2 emissions
          ECO2build <- function(x) {
            x <- mutate(x, ECO2_micro = as.numeric(
              ifelse(hour == 0, abs(-1.57*max(0, 9.9 - airT)-3.6/volume), 
                     ifelse(hour == 1, abs(-1.73*max(0, 12.8 - airT)- 4.2/volume),
                            ifelse(hour == 2, abs(-1.73*max(0,12 - airT)-4.64/volume),
                                   ifelse(hour == 3, abs(-1.5*max(0, 13.1 - airT)-3.44/volume), 
                                          ifelse(hour == 4, abs(-1.38*max(0, 12.4 - airT)-3.86/volume), 
                                                 ifelse(hour == 5, abs(-1.7*max(0, 12.9 - airT)-2.46/volume),
                                                        ifelse(hour == 6, abs(-2.55*max(0, 12.7 - airT)-0.77/volume), 
                                                               ifelse(hour == 7, abs(-2.29*max(0, 10.2 - airT)-2.18/volume), 
                                                                      ifelse(hour == 8, abs(-2.88*max(0, 12.7 - airT)- (-1.5)/volume),
                                                                             ifelse(hour == 9, abs(-2.97*max(0, 12.3 - airT)-(-1.8)/volume), 
                                                                                    ifelse(hour == 10, abs(-1.27*max(0, 13.3 - airT)-0.59/volume), 
                                                                                           ifelse(hour == 11, abs(-1.45*max(0, 15.3 - airT)-0.57/volume),
                                                                                                  ifelse(hour == 12, abs(-1.48*max(0, 15 - airT)-1.1/volume), 
                                                                                                         ifelse(hour == 13, abs(-1.45*max(0, 13.7 - airT)-1.7/volume), 
                                                                                                                ifelse(hour == 14, abs(-1.34*max(0, 13.2 - airT)-2.11/volume),
                                                                                                                       ifelse(hour == 15, abs(-1.35*max(0, 13.5 - airT)-3.22/volume), 
                                                                                                                              ifelse(hour == 16, abs(-1.36*max(0, 13 - airT)-4.39/volume), 
                                                                                                                                     ifelse(hour == 17, abs(-1.29*max(0, 12 - airT)-5.4/volume),
                                                                                                                                            ifelse(hour == 18, abs(-1.44*max(0, 10.6 - airT)-7.3/volume), 
                                                                                                                                                   ifelse(hour == 19, abs(-1.44*max(0, 7 - airT)-9.77/volume), 
                                                                                                                                                          ifelse(hour == 20, abs(-2.7*max(0, 6 - airT)-6.87/volume),
                                                                                                                                                                 ifelse(hour == 21, abs(-1.24*max(0, 8 - airT)-4.29/volume), 
                                                                                                                                                                        ifelse(hour == 22, abs(-0.48*max(0, 10.3 - airT)-3.68/volume),
                                                                                                                                                                               ifelse(hour == 23, abs(-1.69*max(0, 12.4 - airT)-3.6/volume), ""))))))))))))))))))))))))))
            
            return(x)
          }
          
          #Calculate CO2 emissions
          ECO2B <- ECO2build(build_model)
          ECO2T_ras <- ECO2B %>% dplyr::select(x, y, ECO2_micro)
          ECO2T_ras = raster::rasterFromXYZ(xyz = ECO2T_ras,crs = mycrs)
          names(ECO2T_ras) <- paste0("CO2B_", mydate)
          
          return(ECO2T_ras)
        
      }
      
      MapHour <- pbapply::pbapply(ihour, 1, model_hour)
      UHIhour <- unlist(MapHour)
      return(UHIhour)
      
    }
    
    MapDay <- apply(iday, 1, model_day)
    UHIday <- unlist(MapDay)
    return(UHIday)
  }
  
}

job_build <- apply(idates, 1, CO2BuildLCZ) #Apply the function
job_airT_stack <- raster::stack(unlist(job_airT)) #Or get raster stack
qtm(job_airT_stack[[c(19)]]) #plot the map



