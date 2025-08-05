# Working doc for Honours project
# Davina Gifford
# 2025

#load libraries
library(tidyverse)
library(ozmaps)
library(terra)
library(sf)
library(skimr)
library(ncdf4)




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

# map for presentation

ggplot() +
  geom_sf(data = sf_map, alpha = 0.3) +
  geom_sf(data = samp_loc_sf, size = 1, fill = "grey", color = "darkgreen", alpha = 0.5) +
  coord_sf(xlim = c(141, 155), ylim = c(-26, -39)) + # let's zoom in 
  theme_minimal()


# summarise data
skim(samples)


# read mooring data

ncpath <- "D:/HONOURS/DavinaG_2025_Honours/data/"
ncname <- "EAC_filled-daily-distance-depth-gridded-product_20120401-20220727"  
ncfname <- paste(ncpath, ncname, ".nc", sep="")
dname <- "tmp"  # note: tmp means temperature (not temporary)

ncin <- nc_open(ncfname)
print(ncin)



# EAC CCI Anomaly over time -----------------------------------------------



# get anomaly data

anomaly <- tibble(
  trip_month = month_data_t$trip_month,
  eac_cci_clim = month_data_t$eac_cci_clim)


month_data_t$trip_month <- as.Date(month_data_t$trip_month)


# Create the plot with a trend line
ggplot(month_data_t, aes(x = trip_month, y = eac_cci_clim)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "darkred", size = 2) +
  geom_smooth(method = "loess", se = TRUE, color = "black", linetype = "dashed") +
  labs(
    title = "EAC CCI Clim Over Time with Trend Line",
    x = "Trip Month",
    y = "EAC CCI Clim"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


