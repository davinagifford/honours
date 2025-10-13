### scratch.R

### Working doc for Honours project
###
### Created by: Davina Gifford
### Last updated: 2025-08-08

#load libraries
library(tidyverse)
library(ozmaps)
library(terra)
library(sf)
library(skimr)
library(ncdf4)



# Create maps for samples -------------------------------------------------




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



# Validation data ---------------------------------------------------------



# read mooring data

ncpath <- "D:/HONOURS/DavinaG_2025_Honours/data/"
ncname <- "EAC_filled-daily-distance-depth-gridded-product_20120401-20220727"  
ncfname <- paste(ncpath, ncname, ".nc", sep="")
dname <- "tmp"  # note: tmp means temperature (not temporary)

ncin <- nc_open(ncfname)
print(ncin)




# graphs for thesis

plot_predictions(cci_str_model3, condition = "mean_strength", points= 0.25) + 
  labs(
    x = "mean_vel"
  )

plot_predictions(cci_str_model3, condition = "month_clim_str", points= 0.25) +
  labs(
    x = "vel_clim"
  )


plot_predictions(cci_str_model_n3, condition = "mean_strength", points= 0.25) + 
  labs(
  x = "mean_vel"
)

plot_predictions(cci_str_model_n3, condition = "month_clim_str", points= 0.25)+
  labs(
    x = "vel_clim"
  )


plot_predictions(cci_str_model_s3, condition = "mean_strength", points= 0.25) +
  labs(
    x = "mean_vel"
  )

plot_predictions(cci_str_model_s3, condition = "month_clim_str", points= 0.25) +
  labs(
    x = "vel_clim"
  )




# how many unique silk ids in the sample data

samples %>% 
  distinct(silk_id) %>% 
  summarise(n = n())

# how many unique pseg_id in the sample data
samples %>% 
  distinct(pseg_id) %>% 
  summarise(n = n())
19


# look at abundance of species

top_species <- catches %>%
  group_by(species_id) %>%
  summarise(total_abundance = sum(abundance, na.rm = TRUE)) %>%
  arrange(desc(total_abundance)) %>%
  mutate(cum_prop = cumsum(total_abundance) / sum(total_abundance)) %>%
  filter(cum_prop <= 0.89) %>%
  pull(species_id)
top_species
# 78% of total abundance is made up of 15 species