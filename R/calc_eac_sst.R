### calc_eac_sst.R
###
### Calculates EAC monthly nominal SST.
###
### Created: 2023-07-24
### Author: Wayne A. Rochester
### Last updated: 2025-04-28
### Edited by: Davina Gifford

library(lubridate)
library(mgcv)
library(dplyr)
library(readr)

sites <- readRDS(file.path("var", "sst_sites.rds"))
samples_sst <- readRDS(file.path("var", "site_sst.rds"))

samples_sst <-
    samples_sst %>%
    rename(sample_time = sst_date) %>%
    left_join(sites, by = "site_id")

sample_stats <-
    samples_sst %>%
    summarise(time0 = floor_date(min(sample_time), unit = "day"),
              mean_time = mean(sample_time),
              mean_lat = mean(latitude))

time0 <- pull(sample_stats, time0)
mean_lat <- pull(sample_stats, mean_lat)
mean_time <- pull(sample_stats, mean_time)

samples_sst <-
    samples_sst %>%
    mutate(doy = yday(sample_time),
           time_x = time_length(interval(time0, sample_time), unit = "day"))

lm_fit <- gam(sst ~ s(doy, bs = "cc") + s(latitude) + s(time_x),
              knots=list(doy = c(0, 365)),
              data = samples)

print(summary(lm_fit))

terms_pred <- predict(lm_fit, type = "terms", newdata = samples_sst)

samples_sst <-
    samples_sst %>%
    mutate(lat_eff = terms_pred[, "s(latitude)"],
           eac_sst = sst - lat_eff) %>%
    select(!lat_eff)

lm_fit2 <- gam(eac_sst ~ s(doy, bs = "cc") + s(latitude) + s(time_x),
               knots = list(doy=c(0, 365)),
               data = samples_sst)
print(summary(lm_fit2))

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
           eac_sst = doy_eff + intercept) %>%
    select(!c(doy_eff, intercept))

month_data <-
    samples_sst %>%
    mutate(sample_month = floor_date(sample_time, unit = "month")) %>%
    group_by(sample_month) %>%
    summarise(sst = mean(sst),
              eac_sst = mean(eac_sst),
              num_samples = n(),
              .groups = "drop")

cat("\nEffect of SST adjustment on monthly averages (should be near zero):\n")
print(month_data %>% reframe(diff_range = range(eac_sst - sst)))

sst_anom_data <- list(lm_fit = lm_fit,
                      samples = samples,
                      climatology = climatology,
                      month_data = month_data)
saveRDS(sst_anom_data, file.path("var", "eac_sst.rds"))
