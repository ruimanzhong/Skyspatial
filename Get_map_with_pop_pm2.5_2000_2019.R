library(readxl)
library(sf)
library(ggplot2)
library(rgeoboundaries)
library(rnaturalearth)
library(SpatialEpi)
library(spdep)
library(INLA)
library(sp)
library(Matrix)
library(spData)
library(akima)
library(terra)
library(raster)
library(leaflet)
library(shiny)
library(rmarkdown)
library(htmltools)
library(tidyr)
library(dplyr)
library(webshot)
library(htmlwidgets)
library(pander)
library(here)

rm(list = ls())
theme_set(theme_minimal())

map_list <- list()

# Get the map of Saudi Arabia
initial_map <- ne_states(country = "Saudi Arabia", returnclass = "sf") %>% 
  mutate(name = replace(name, name == name[5], "Asir")) %>% 
  st_set_crs(4326) %>% 
  arrange(name) %>% 
  dplyr::select(name, geometry)
bb <- st_bbox(initial_map)
  
  # Loop through the years 2000 to 2019
  for (year in 2000:2019) {
    
    # Load and process population raster
    pop_file <- here(paste0("SApop/sau_pd_", year, "_1km.tif"))
    r <- terra::rast(pop_file) %>% 
      terra::project("EPSG:4326") %>% 
      terra::aggregate(fact = 10, fun = sum)
    
    # Load and process PM2.5 raster
    pm_file <- here(paste0("SApm2.5/pm2.5.", year, ".tif"))
    data <- terra::rast(pm_file) %>% 
      terra::aggregate(fact = 10, fun = mean) %>%
      crop(bb) %>%
      terra::project("EPSG:4326")
    
    # Update map with extracted values for population and PM2.5
    map <- initial_map %>%
      mutate(
        population = terra::extract(r, ., fun = sum, na.rm = TRUE)[[2]],
        pm25 = terra::extract(data, ., fun = mean, na.rm = TRUE)[[2]]
      )
    
    # Convert PM2.5 raster to a spatial data frame
    dpol <- data %>%
      as.data.frame(xy = TRUE) %>%
      rename(lon = x, lat = y, value = paste0("pm2.5.", year)) %>%
      st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
      st_filter(map) %>%
      mutate(pm25g = as.numeric(value),
             pm25Xpop = pm25g * terra::extract(r, st_coordinates(.))[[1]],
             pop = terra::extract(r, st_coordinates(.))[[1]]
      )
    
    # Spatial intersection between map and PM2.5 spatial data
    inter <- st_intersects(map, dpol)
    
    # Calculate weighted population PM2.5 values
    map$pm25xpop <- sapply(inter, FUN = function(x){sum(dpol[x, ]$pm25Xpop, na.rm = TRUE)})
    map$sumpop <- sapply(inter, FUN = function(x){sum(dpol[x, ]$pop, na.rm = TRUE)})
    map$WPM2.5 <- map$pm25xpop / map$sumpop
    
    # Select relevant columns and add the current year
    map <- map %>%
      dplyr::select(name, population, WPM2.5) %>%
      mutate(year = year)
    
    # Store the map for this year in the list
    map_list[[as.character(year)]] <- map
  }

# Combine all yearly maps into one data frame
map <- bind_rows(map_list)

# The following is just to get clipped.raster.dpol

coord.dpol <- st_coordinates(dpol)

raster.dpol <- data.frame(lon = coord.dpol[,1], lat = coord.dpol[,2], value = dpol$pm25g) %>% 
  with(interp(x = lon, y = lat, z = value, xo = seq(min(lon), max(lon), length = 100), yo = seq(min(lat), max(lat), length = 100))) %>%
  raster()

# Load the Saudi Arabia boundary
boundaryregion_sp <- geoboundaries("Saudi Arabia") %>% st_set_crs(4326) %>% st_transform(crs(proj4string(raster.dpol))) %>% as("Spatial")

# Clip the raster to the Saudi Arabia boundary
clipped.raster.dpol <- mask(raster.dpol, boundaryregion_sp)


# Save the final map as an RDS file
save(map, clipped.raster.dpol, file = here("map_pop_pm2.5_2000_2019.RData"))
