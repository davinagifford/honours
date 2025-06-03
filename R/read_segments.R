### read_segments.R
###
### Reads and selects CPR segments for use in the analysis.
###
### Created: 2023-07-19
### Author: Wayne A. Rochester

library(sf)
library(dplyr)
library(readr)

segments <- read_csv(file.path("data", "samples.csv"),
                     col_types = "iiiddTTicicccicicd")

saveRDS(segments, file.path("var", "segments_all.rds"))

mask_lyr <- read_sf(file.path("shape", "eac_mask.shp"), crs = 4326)
mask_lyr <- mask_lyr %>% select(cat)

segs_t <-
    segments %>%
    select(pseg_id, latitude, longitude) %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
    st_join(mask_lyr) %>%
    st_drop_geometry() %>%
    filter(!is.na(cat))

segments <-
    segments %>%
    semi_join(segs_t, by = "pseg_id") %>%
    filter(latitude > -37.67) # NSW / Vic border

saveRDS(segments, file.path("var", "segments.rds"))
