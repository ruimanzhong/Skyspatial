


library(readxl)
library(sf)
library(ggplot2)
library(rgeoboundaries)
library(rnaturalearth)
library(dplyr)
library(SpatialEpi)
library(spdep)
library(INLA)
library(terra)
library(tigris)



##################
# Get the map of Saudi Arabia
##################

map <- ne_states(country = "Saudi Arabia", returnclass = "sf")
map$name[5]="Asir" #change spellings
map <- st_set_crs(map, 4326) #certain projection
map <- map[order(map$name),] #alphabetically ordered
map<-map %>% dplyr::select(name,geometry)
maplatlong <- st_transform(map, 4326)




#################################
# Population of Saudi Arabia
#################################

# Get the population distribution of Saudi Arabia in year of 2017 from here
# https://hub.worldpop.org/geodata/listing?id=76 



# # Read the saudi arabia population raster for year 2018
r <- terra::rast("~/Desktop/files/PhD/Projects/project3/SummerProject3/Nojoud/Data2017/data/pop2017.tif")
# # Ensure raster is in EPSG:4326 
r <- terra::project(r, "EPSG:4326")
plot(r)
# 
# # Transform long/lat to UTM (same as map)
r <- terra::project(r, "EPSG:4326") 
# 
# # # Sum population in regions of map
map$popraster <- terra::extract(r, map, sum, na.rm = TRUE)$pop2017
# # 
# # # modify the data for plot
r2 <- as.data.frame(r, xy = TRUE)
names(r2) <- c("lon", "lat", "value") 
r2 <- st_as_sf(r2, coords = c("lon", "lat"), crs = 4326)




# read the pm2.5 data as grid for 2017

data <- terra::rast("~/Desktop/files/PhD/Projects/project3/SummerProject3/Nojoud/Data2017/data/pm2.5.2017.tif")
data <- terra::aggregate(data, fact=2)

bb <- st_bbox(maplatlong)
data <- crop(data, bb)
plot(data)
# Transform long/lat to UTM (same as map)
data <- terra::project(data, "EPSG:4326") 

# Sum pm in regions of map
map$pm2.5 <- terra::extract(data, map, mean, na.rm = TRUE)$pm2.5.2017

# modify the data for plot
data2 <- as.data.frame(data, xy = TRUE)
names(data2) <- c("lon", "lat", "value") 


dpol <- st_as_sf(data2, coords = c("lon", "lat"), crs = 4326)
dpol <- st_filter(dpol, map)
dpol$pm252012g <- as.numeric(dpol$value)






##################
# weighted PM25 
##################


# PM25 weighted by population pm25 * pop / sum(pop)
# (PM25 at the coordinates PM25 raster) * (Population at coordinates PM25 raster)

dpol$pm25Xpop <- as.vector(dpol$pm252012g * terra::extract(r, st_coordinates(dpol)))[[1]]
# (Population at coordinates PM25 raster)

dpol$pop <- as.vector(terra::extract(r, st_coordinates(dpol)))[[1]]
# Mean air pollution within counties weighted by population
inter <- st_intersects(map, dpol)
map$pm25xpop <- sapply(inter, FUN = function(x){sum(dpol[x, ]$pm25Xpop, na.rm = TRUE)})
map$sumpop <- sapply(inter, FUN = function(x){sum(dpol[x, ]$pop, na.rm = TRUE)})
map$pm25weightedpop <- map$pm25xpop/map$sumpop
