### read_catch_data.R
###
### Reads the catch data.
###
### Created: 2023-07-19
### Author: Wayne A. Rochester

library(lubridate)
library(tidyverse)
library(dplyr)
library(tidyr)
library(readr)

segments <- readRDS(file.path("var", "segments.rds"))

samples <- read_csv(file.path("data", "samples.csv"),
                    col_types = "iiiddTTicicccicicd")
species <- read_csv(file.path("data", "species.csv"), col_types = "ic")
catches <- read_csv(file.path("data", "catches.csv"), col_types = "iid")

## Taxononic identification was changing significantly before April 2011.

samples <-
    samples %>%
    semi_join(segments, by = "pseg_id") %>%
    filter(process_date >= ymd_hms("2011-04-01 00:00:00",
                                   tz = tz(process_date))) %>%
    semi_join(catches %>% distinct(pseg_id), by = "pseg_id")
catches <-
    catches %>%
    semi_join(samples, by = "pseg_id")
species <-
    species %>%
    semi_join(catches %>% distinct(species_id), by = "species_id")

catch_data <- list(samples = samples,
                   species = species,
                   catches = catches)

saveRDS(catch_data, file.path("var", "catch_data.rds"))
