### calc_eac_cci.R
###
### Calculates the EAC copepod composition index.
###
### Created: 2023-07-25
### Author: Wayne A. Rochester
### Last updated: 2025-04-28
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
library(gratia)


catch_data <- readRDS(file.path("var", "catch_data.rds"))

samples <- catch_data$samples
species <- catch_data$species
catches <- catch_data$catches

catch_mtx <- catches %>%
  mutate(
    pseg_id = factor(pseg_id, levels = samples$pseg_id), # create factors of the sample ID
    species_id = factor(species_id, levels = species$species_id)
  ) %>% # create factors of the species ID
  group_by(pseg_id, species_id) %>% # grouping by sample, then by species
  summarise(max_abundance = max(abundance, na.rm = TRUE), .groups = "drop") %>% # creating a column for max abundance
  pivot_wider( # create the matrix
    names_from = species_id,
    values_from = max_abundance,
    values_fill = 0,
    id_expand = TRUE, # ensures that all combinations of pseg_id and species_id are included
    names_expand = TRUE, 
  )


# Convert to matrix, excluding the grouping column
catch_mtx_numeric <- as.matrix(catch_mtx[,-1])
# make row names the sample id 
rownames(catch_mtx_numeric) <- catch_mtx$pseg_id

stopifnot(rownames(catch_mtx_numeric) == as.character(samples$pseg_id))
stopifnot(colnames(catch_mtx_numeric) == as.character(species$species_id))

# double square root transformation followed by Hellinger transformation
catch_mtx_trfm <- sqrt(sqrt(catch_mtx_numeric))
catch_mtx_trfm <- decostand(catch_mtx_trfm, "hellinger")



# RDA ---------------------------------------------------------------------

# RDA on the transformed data against 

rda_fit <- rda(catch_mtx_trfm ~ sst, data=samples) 

summary(rda_fit)

print(summary(rda_fit)$cont$importance[, 1:4])


# Test significance

anova.cca(rda_fit, step = 1000)
RsquareAdj(rda_fit)


# Extract scores

samp_score <- scores(rda_fit, display="wa")[, 1] # Sample (site) scores on RDA1. wa = weighted average
sp_score <- scores(rda_fit, display="sp")[, 1] # Species scores on RDA1

# convert scores to positive, for interpretability
# ensure higher score = higher SST

if (cor(samples$sst, samp_score, use="complete.obs") < 0) { # checks correlation between SST and Sample score
    samp_score <- -samp_score
    sp_score <- -sp_score
}

# add scores to the dataframes, for use in analysis

samples <- samples %>% mutate(rda_score = samp_score)
#species <- species %>% mutate(rda_score = sp_score)

# ensure that species scores have same order as the species table
sp_score_df <-
  sp_score %>%
  enframe(name = "species_id", value = "rda_score") %>%
  mutate(species_id = as.integer(species_id))

species <- species %>% left_join(sp_score_df, by = "species_id")


# summarise stats of samples 
sample_stats <-
    samples %>%
    summarise(time0 = floor_date(min(sample_time), unit = "day"), # starting date
              mean_time = mean(sample_time), # average sample time
              mean_lat = mean(latitude)) # average sample latitude - create a fixed-location "typical sample" in the climatology model.

time0 <- pull(sample_stats, time0)
mean_lat <- pull(sample_stats, mean_lat)
mean_time <- pull(sample_stats, mean_time)

samples <-
    samples %>%
    mutate(doy = yday(sample_time), # define Day of Year
           time_x = time_length(interval(time0, sample_time), unit = "day")) # provides a numeric variable suitable for smoothing in a GAM 



# GAM ---------------------------------------------------------------------



lm_fit <- gam(rda_score ~ s(doy, bs = "cc", k = 5) + s(latitude) +
                  s(time_x, k = 25),
              knots = list(doy = c(0, 365)), # sets the bounds on doy - limit to 365 days
              data = samples)

#check the gam
gam.check(lm_fit)


draw(lm_fit)



print(summary(lm_fit))



# returns the contribution of each individual smooth term for each row.
# Enables removal of confounding terms (like latitude) from your response

terms_pred <- predict(lm_fit, type = "terms", newdata = samples)
head(terms_pred)



samples <-
    samples %>%
    mutate(lat_eff = terms_pred[, "s(latitude)"],  # extract the latitude effect 
           eac_cci = rda_score - lat_eff) %>% # calculate CCI by removing latitude effect from RDA Score
    select(!lat_eff) # remove latitude effect column 


## Calculate a daily climatology of the EAC CCI as the sum of the
## relevant GAM prediction components (the intercept and day of year).

climatology <-
    tibble(sample_time = floor_date(mean_time, unit = "year") + days(0:364),
           latitude = mean_lat,
           doy = yday(sample_time),
           time_x = time_length(interval(time0, mean_time), unit = "day"))

clim_terms_pred <- predict(lm_fit, type = "terms", newdata = climatology)

head(clim_terms_pred)

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
