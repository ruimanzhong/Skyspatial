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


rm(list = ls())
theme_set(theme_minimal())

# Get the map of Saudi Arabia
map <- ne_states(country = "Saudi Arabia", returnclass = "sf") %>% 
  mutate(name = replace(name, name == name[5], "Asir")) %>% 
  st_set_crs(4326) %>% 
  arrange(name) %>% 
  dplyr::select(name, geometry)
bb <- st_bbox(map)

r <- terra::rast(here("SApop/sau_pd_2000_1km.tif")) %>% 
  terra::project("EPSG:4326") %>% 
  terra::aggregate(fact = 10, fun = sum)

data <- terra::rast(here("SApm2.5/pm2.5.2000.tif")) %>% 
  terra::aggregate(fact = 10, fun = mean) %>%
  crop(bb) %>%
  terra::project("EPSG:4326")

map <- map %>%
  mutate(
    popraster = terra::extract(r, ., fun = sum, na.rm = TRUE)$sau_pd_2000_1km,
    pm2.5 = terra::extract(data, ., fun = mean, na.rm = TRUE)$pm2.5.2000
  )



dpol <- data %>%
  as.data.frame(xy = TRUE) %>%
  rename(lon = x, lat = y, value = pm2.5.2000) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_filter(map) %>%
  mutate(pm25g = as.numeric(value)) %>%
  mutate(
    pm25Xpop = pm25g * terra::extract(r, st_coordinates(.))[[1]],
    pop = terra::extract(r, st_coordinates(.))[[1]]
  )

inter <- st_intersects(map, dpol)

map$pm25xpop <- sapply(inter, FUN = function(x){sum(dpol[x, ]$pm25Xpop, na.rm = TRUE)})
map$sumpop <- sapply(inter, FUN = function(x){sum(dpol[x, ]$pop, na.rm = TRUE)})
map$pm25weightedpop <- map$pm25xpop/map$sumpop

map <- map %>%
  dplyr::select(name, popraster, pm25weightedpop) %>%
  mutate(year = 2000)
