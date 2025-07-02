### calc_eac_cci_south.R
###
### Calculates the EAC copepod composition index for the EAC below the separation zone at 32 degrees south
###
### Created: 2023-07-25
### Author: Wayne A. Rochester
### Last updated: 2025-07-02
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
samples_south <- samples %>% 
  filter(latitude >= "-32")

# create the grouping index for organising the data
# tapply_index <-
#    catches %>%
#    mutate(pseg_id = factor(pseg_id, levels = samples_north$pseg_id), # create factors of the sample ID
#           species_id = factor(species_id, levels = species$species_id)) %>% # create factors of the species ID
#    select(pseg_id, species_id)

# Tapply - base R
# Find the maximum abundance for each combination of sample and species and creates a matrix

#catch_mtx <- tapply(catches$abundance, tapply_index, max, default = 0) 

# do the same but use tidyverse code


catch_mtx <- catches %>%
  mutate(
    pseg_id = factor(pseg_id, levels = samples_south$pseg_id), # create factors of the sample ID
    species_id = factor(species_id, levels = species$species_id)
  ) %>% # create factors of the species ID
  group_by(pseg_id, species_id) %>% # grouping by sample, then by species
  summarise(max_abundance = max(abundance, na.rm = TRUE), .groups = "drop") %>% # creating a column for max abundance
  pivot_wider( # create the matrix
    names_from = species_id,
    values_from = max_abundance,
    values_fill = 0
  )


# Convert to matrix, excluding the grouping column
catch_mtx_numeric <- as.matrix(catch_mtx[,-1])
# make row names the sample id 
rownames(catch_mtx_numeric) <- catch_mtx$pseg_id





# double square root transformation followed by Hellinger transformation
catch_mtx_trfm <- sqrt(sqrt(catch_mtx_numeric))
catch_mtx_trfm <- decostand(catch_mtx_trfm, "hellinger")


# Ensure rownames match the sample IDs
rownames(catch_mtx_trfm) <- rownames(catch_mtx_numeric) <- catch_mtx$pseg_id


# Remove rows with NA rownames
catch_mtx_trfm <- catch_mtx_trfm[!is.na(rownames(catch_mtx_trfm)), ]


# Subset samples_south to match the matrix rows
samples_south_aligned <- samples_south %>%
  filter(pseg_id %in% rownames(catch_mtx_trfm)) %>%
  arrange(factor(pseg_id, levels = rownames(catch_mtx_trfm)))



# RDA using the transformed catch data and the SST values from the sample df
rda_fit_south <- rda(catch_mtx_trfm ~ sst, data=samples_south_aligned)

print(summary(rda_fit_south)$cont$importance[, 1:4])

# check significance of RDA model
anova_rda <- anova.cca(rda_fit_south, permutations = 999)
print(anova_rda)

samp_score_s <- scores(rda_fit_south, display="wa")[, 1]
sp_score_s <- scores(rda_fit_south, display="sp")[, 1]

if (cor(samples_south$sst, samp_score_s, use="complete.obs") < 0) {
    samp_score_s <- -samp_score_s
    sp_score_s <- -sp_score_s
}

samples_south <- samples_south %>% mutate(rda_score = samp_score_s)

head(samples_south)


species_s <- species %>% mutate(rda_score = sp_score_s)


# get summary stats
sample_stats_s <-
    samples_south %>%
    summarise(time0 = floor_date(min(sample_time), unit = "day"),
              mean_time = mean(sample_time),
              mean_lat = mean(latitude))

time0 <- pull(sample_stats_s, time0)
mean_lat <- pull(sample_stats_s, mean_lat)
mean_time <- pull(sample_stats_s, mean_time)

# add in temporal data to the samples_south df
samples_south <-
    samples_south %>%
    mutate(doy = yday(sample_time),
           time_x = time_length(interval(time0, sample_time), unit = "day"))

# GAM
# cyclic cubic splines - these make the start and end f year connect smoothly
# k is the number of 'knots' to make the spline. Small ,5 = simple, less wiggly curve
lm_fit_south <- gam(rda_score ~ s(doy, bs = "cc", k = 5) + # day of year, , K?
                s(latitude) + # latitude smoother
                s(time_x, k = 8), # time
              knots = list(doy = c(0, 365)), # boundary for the day cycle
              data = samples_south)

print(summary(lm_fit_south))

terms_pred <- predict(lm_fit_south, type = "terms", newdata = samples_south)

samples_south <-
    samples_south %>%
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

clim_terms_pred <- predict(lm_fit_south, type = "terms", newdata = climatology)

climatology <-
    climatology %>%
    mutate(doy_eff = clim_terms_pred[, "s(doy)"],
           intercept = coef(lm_fit_south)["(Intercept)"],
           eac_cci = doy_eff + intercept) %>%
    select(!c(doy_eff, intercept))


## When calculating monthly averages, we use average trip times rather
## than sample times to ensure that all samples from one trip are
## assigned to the same month. (Trips are only a few days long.)

south_data <-
    samples_south %>%
    group_by(trip_id) %>%
    mutate(trip_time = mean(sample_time)) %>%
    ungroup() %>%
    mutate(trip_month = floor_date(trip_time, unit = "month")) %>%
    group_by(trip_month) %>%
    summarise(eac_cci = mean(eac_cci),
              num_samples = n(),
              .groups = "drop")

cci_data_south <- list(samples = samples_south,
                 species = species_s,
                 catches = catches,
                 catch_mtx = catch_mtx,
                 catch_mtx_trfm = catch_mtx_trfm,
                 rda_fit = rda_fit_south,
                 lm_fit = lm_fit_south,
                 climatology = climatology,
                 month_data = south_data)
saveRDS(cci_data_south, file.path("var", "eac_cci_south.rds"))
