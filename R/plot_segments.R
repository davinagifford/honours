### plot_segments.R
###
### Plots the CPR segments.
###
### Created: 2023-07-19
### Author: Wayne A. Rochester

library(sf)
library(dplyr)
library(ggplot2)

segments <- readRDS(file.path("var", "segments.rds"))
segments_all <- readRDS(file.path("var", "segments_all.rds"))
coast_lyr <- read_sf(file.path("shape", "aus_region.shp"), crs = 4326)
mask_lyr <- read_sf(file.path("shape", "eac_mask.shp"), crs = 4326)

xlim <- c(140, 160)
ylim <- c(-45, -25)
xlim <- xlim + c(-0.5, 0.5)
ylim <- ylim + c(-0.5, 0.5)

p <-
    ggplot() +
    geom_sf(data = coast_lyr, colour = "black", fill = "gray") +
    geom_sf(data = mask_lyr, colour = "blue", fill = "cyan") +
    geom_point(data = segments_all, mapping = aes(longitude, latitude),
               colour = "gray") +
    geom_point(data = segments, mapping = aes(longitude, latitude),
               colour = "red") +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(x = NULL, y = NULL)
ggsave(file.path("output", "segment_map.png"),
       plot = p, width = 680 / 96, height = 680 / 96, dpi = 96,
       device = png)

p <-
    ggplot(segments, aes(sample_time, latitude, colour = route)) +
    geom_point()
ggsave(file.path("output", "segment_ts.png"),
       plot = p, width = 800 / 96, height = 450 / 96, dpi = 96,
       device = png)

