library(here)
library(terra)
library(dplyr)
library(sf)
rm(list = ls())

# Get the map of Saudi Arabia
map <- ne_states(country = "Saudi Arabia", returnclass = "sf") |> 
  mutate(name = replace(name, name == name[5], "Asir")) |> 
  st_set_crs(4326) |> arrange(name) |> 
  dplyr::select(name, geometry)

# List of raster files
raster_files <- list.files(path = "SApop", pattern = "sau_pd_.*_1km.tif", full.names = TRUE)

# Initialize a column to store population data for each year
for (raster_file in raster_files) {
  # Extract the year from the filename
  year <- sub("sau_pd_(\\d{4})_1km.tif", "\\1", basename(raster_file))
  
  # Load and process the raster file
  r <- rast(raster_file) |> 
    project("EPSG:4326") |> 
    aggregate(fact = 10, fun = sum)
  
  # Extract the variable name from the raster
  var_name <- names(r)
  
  # Extract population data
  pop_data <- terra::extract(r, map, sum, na.rm = TRUE)
  
  # Add population data to the map data frame
  map[[paste0("pop_", year)]] <- pop_data[[var_name]]
}

# Check the updated map data frame
print(map)

