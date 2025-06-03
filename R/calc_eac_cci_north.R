### calc_eac_cci_north.R
###
### Calculates the EAC copepod composition index.
###
### Created: 2023-07-25
### Author: Wayne A. Rochester
### Last updated: 2025-05-12
### Edited by: Davina Gifford

## The EAC copepod composition index (EAC CCI) is calculated by an RDA
## followed by a GAM. The RDA is of copepod composition on SST. The
## EAC CCI is calculated as the RDA sample score (first axis) minus
## the latitude effect. The latitude effect is estimated with the GAM.

library(lubridate)
library(dplyr)
library(tidyr)
library(mgcv)
library(vegan)
library(tidyverse)


catch_data <- readRDS(file.path("var", "catch_data.rds"))

samples <- catch_data$samples
species <- catch_data$species
catches <- catch_data$catches

# filter to north of separation
samples_north <- samples %>% 
  filter(latitude <= "-31")

# create the grouping index for organising the data
tapply_index <-
    catches %>%
    mutate(pseg_id = factor(pseg_id, levels = samples_north$pseg_id), # create factors of the sample ID
           species_id = factor(species_id, levels = species$species_id)) %>% # create factors of the species ID
    select(pseg_id, species_id)

# Tapply - base R
# Find the maximum abundance for each combination of sample and species and creates a matrix

catch_mtx <- tapply(catches$abundance, tapply_index, max, default = 0) 

# double square root transformation followed by Hellinger transformation
catch_mtx_trfm <- sqrt(sqrt(catch_mtx))
catch_mtx_trfm <- decostand(catch_mtx_trfm, "hellinger")

# RDA using the transformed catch data and the SST values from the sample df
rda_fit_north <- rda(catch_mtx_trfm ~ sst, data=samples_north)

print(summary(rda_fit_north)$cont$importance[, 1:4])


samp_score_n <- scores(rda_fit_north, display="wa")[, 1]
sp_score_n <- scores(rda_fit_north, display="sp")[, 1]

if (cor(samples_north$sst, samp_score_n, use="complete.obs") < 0) {
    samp_score_n <- -samp_score_n
    sp_score_n <- -sp_score_n
}

samples_north <- samples_north %>% mutate(rda_score = samp_score_n)

head(samples_north)


species_n <- species %>% mutate(rda_score = sp_score_n)


# get summary stats
sample_stats_n <-
    samples_north %>%
    summarise(time0 = floor_date(min(sample_time), unit = "day"),
              mean_time = mean(sample_time),
              mean_lat = mean(latitude))

time0 <- pull(sample_stats_n, time0)
mean_lat <- pull(sample_stats_n, mean_lat)
mean_time <- pull(sample_stats_n, mean_time)

# add in temporal data to the samples_north df
samples_north <-
    samples_north %>%
    mutate(doy = yday(sample_time),
           time_x = time_length(interval(time0, sample_time), unit = "day"))

# GAM
# cyclic cubic splines - these make the start and end f year connect smoothly
# k is the number of 'knots' to make the spline. Small ,5 = simple, less wiggly curve
lm_fit <- gam(rda_score ~ s(doy, bs = "cc", k = 5) + # day of year, , K?
                s(latitude) + # latitude smoother
                s(time_x, k = 8), # time
              knots = list(doy = c(0, 365)), # boundary for the day cycle
              data = samples_north)

print(summary(lm_fit))

terms_pred <- predict(lm_fit, type = "terms", newdata = samples)

samples <-
    samples %>%
    mutate(lat_eff = terms_pred[, "s(latitude)"],
           eac_cci = rda_score - lat_eff) %>%
    select(!lat_eff)


## Calculate a daily climatology of the EAC CCI as the sum of the
## relevant GAM prediction components (the intercept and day of year).

climatology <-
    tibble(sample_time = floor_date(mean_time, unit = "year") + days(0:364),
           latitude = mean_lat,
           doy = yday(sample_time),
           time_x = time_length(interval(time0, mean_time), unit = "day"))

clim_terms_pred <- predict(lm_fit, type = "terms", newdata = climatology)

climatology <-
    climatology %>%
    mutate(doy_eff = clim_terms_pred[, "s(doy)"],
           intercept = coef(lm_fit)["(Intercept)"],
           eac_cci = doy_eff + intercept) %>%
    select(!c(doy_eff, intercept))


## When calculating monthly averages, we use average trip times rather
## than sample times to ensure that all samples from one trip are
## assigned to the same month. (Trips are only a few days long.)

month_data <-
    samples %>%
    group_by(trip_id) %>%
    mutate(trip_time = mean(sample_time)) %>%
    ungroup() %>%
    mutate(trip_month = floor_date(trip_time, unit = "month")) %>%
    group_by(trip_month) %>%
    summarise(eac_cci = mean(eac_cci),
              num_samples = n(),
              .groups = "drop")

cci_data <- list(samples = samples,
                 species = species,
                 catches = catches,
                 catch_mtx = catch_mtx,
                 catch_mtx_trfm = catch_mtx_trfm,
                 rda_fit = rda_fit,
                 lm_fit = lm_fit,
                 climatology = climatology,
                 month_data = month_data)
saveRDS(cci_data, file.path("var", "eac_cci.rds"))
