# Working doc for Honours project
# Davina Gifford
# 2025

#load libraries
library(tidyverse)
library(ozmaps)
library(terra)
library(sf)
library(skimr)




samples %>% 
  distinct(trip_id) %>% 
  summarise(n = n())

samples %>% 
  distinct(silk_id) %>% 
  summarise(n = n())


# map the samples

sample_location <- samples %>% 
  select(latitude:longitude)

sf_map <- ozmap("states")

# make the sample locations a simple feature
samp_loc_sf <- st_as_sf(sample_location, coords = c("longitude", "latitude"))
st_crs(samp_loc_sf) <-"EPSG:4283"
st_crs(sf_map)

# make the locations have the same crs as the map
samp_loc_sf <- st_transform(samp_loc_sf, st_crs(sf_map))


# make a map

ggplot() +
  geom_sf(data = sf_map, alpha = 0.3) +
  geom_sf(data = samp_loc_sf, size = 1, fill = "darkgreen") +
  coord_sf(xlim = c(135, 160), ylim = c(-10, -45)) + # let's zoom in 
  theme_minimal()


# summarise data
skim(samples)