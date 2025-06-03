### plot_sst_sites.R
###
### Plots a map of the SST sampling sites.
###
### Created: 2023-07-24
### Author: Wayne A. Rochester

library(sf)
library(dplyr)
library(ggplot2)

xlim <- c(139.8, 157.5)
ylim <- c(-39.9, -25.1)

sst_sites <- readRDS(file.path("var", "sst_sites.rds"))
segments <- readRDS(file.path("var", "segments.rds"))
coast_lyr <- read_sf(file.path("shape", "aus_region.shp"), crs = 4326)

p <-
    ggplot() +
    geom_sf(data = coast_lyr, colour = "black", fill = "gray") +
    geom_point(data = segments, mapping = aes(longitude, latitude),
               colour = "gray") +
    geom_point(data = sst_sites, mapping = aes(longitude, latitude),
               colour = "red") +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(x = NULL, y = NULL)
ggsave(file.path("output", "sst_sites.png"),
       plot = p, width = 680 / 96, height = 680 / 96, dpi = 96,
       device = png)

