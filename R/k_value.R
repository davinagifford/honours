# choosing k values for the GAM
#
# including validation and visualisations
# 
# 19 May 2025
#
# Davina Gifford



# load libraries

library(mgcv)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

# Define k values to test
k_grid <- expand.grid(
  k_doy = c(5, 10, 15, 20, 25),
  k_time = c(5, 10, 15, 20, 25)
)

# Fit models and extract metrics
model_results <- k_grid %>%
  mutate(
    model = map2(k_doy, k_time, ~ gam(
      rda_score ~ s(doy, bs = "cc", k = .x) +
        s(latitude, k = 9) +
        s(time_x, k = .y),
      knots = list(doy = c(0, 365)),
      data = samples
    )),
    gcv = map_dbl(model, ~ .x$gcv.ubre),
    aic = map_dbl(model, AIC),
    edf_doy = map_dbl(model, ~ summary(.x)$s.table["s(doy)", "edf"]),
    edf_time = map_dbl(model, ~ summary(.x)$s.table["s(time_x)", "edf"])
  )

ggplot(model_results, aes(k_doy, k_time, fill = gcv)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c() +
  labs(title = "GCV score across k values", x = "k for doy", y = "k for time_x")

ggplot(model_results, aes(k_doy, k_time, fill = aic)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c() +
  labs(title = "AIC across k values", x = "k for doy", y = "k for time_x")

model_results %>%
  arrange(gcv) %>%
  slice(1:3)
