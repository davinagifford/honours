### plot_eac_cci.R
###
### Plots the EAC copepod composition index.
###
### Created: 2023-07-21
### Author: Wayne A. Rochester
### Last updated: 2025-04-28
### Edited by: Davina Gifford

library(lubridate)
library(vegan)
library(mgcv)
library(dplyr)
library(tidyr)
library(ggplot2)

calc_seas_vars <- function (x) {
    ## Calculates season as x-y components for use as regression variables.

    doy <- pmin(yday(x), 365)
    orig_doy <- yday(ymd_hms("2019-01-15 00:00:00"))

    theta <- -(doy - orig_doy) / 365 * pi * 2

    tibble(win_sum = cos(theta),
           aut_spr = sin(theta))
}


cci_data <- readRDS(file.path("var", "eac_cci.rds"))

samples <- cci_data$samples
rda_fit <- cci_data$rda_fit
lm_fit <- cci_data$lm_fit
climatology <- cci_data$climatology
month_data <- cci_data$month_data


month_data_rf <- samples %>%
  group_by(trip_id) %>%
  mutate(trip_time = mean(sample_time)) %>%
  ungroup() %>%
  mutate(trip_month = floor_date(trip_time, unit = "month")) %>%
  group_by(trip_month) %>%
  summarise(eac_cci_rf = mean(eac_cci_rf_pca_latadj, na.rm = TRUE),
            num_samples = n(), .groups = "drop")



month_data_t_rf <- month_data_rf %>%
  mutate(sample_time = trip_month + ddays(days_in_month(trip_month)) * 0.5,
         sample_time = floor_date(sample_time, unit = "day"),
         doy = yday(sample_time)) %>%
  left_join(climatology %>% select(doy, eac_cci_clim = eac_cci), by = "doy") %>%
  mutate(anom_label = if_else(eac_cci_rf >= eac_cci_clim, "Anomaly (+)", "Anomaly (-)"))

time_range <- range(month_data_t_rf$sample_time)

clim_points_rf <- tibble(sample_time = seq(time_range[1], time_range[2], by = "1 day")) %>%
  mutate(doy = pmin(yday(sample_time), 365)) %>%
  left_join(climatology %>% select(doy, eac_cci), by = "doy")



p_rf <-
  ggplot(mapping = aes(sample_time, eac_cci_rf)) +
  geom_segment(data = month_data_t_rf,
               mapping = aes(xend = sample_time, yend = eac_cci_clim,
                             colour = anom_label, linetype = anom_label)) +
  geom_line(data = clim_points_rf,
            mapping = aes(y = eac_cci,
                          linetype = "Climatology",
                          colour = "Climatology")) +
  geom_point(data = month_data_t_rf, size = 2, shape = 21, fill = "gray") +
  scale_colour_manual(breaks = c("Anomaly (+)", "Anomaly (-)", "Climatology"),
                      values = c("red", "blue", "black")) +
  scale_linetype_manual(breaks = c("Anomaly (+)", "Anomaly (-)", "Climatology"),
                        values = c("solid", "dashed", "dotted")) +
  scale_x_datetime(date_breaks = "1 year", minor_breaks = NULL,
                   date_labels = "%Y") +
  labs(x = "Time",
       y = "EAC copepod composition index (RF + PCA)",
       colour = NULL,
       linetype = NULL,
       title = "EAC copepod composition index (RF + PCA, monthly average)")

ggsave(file.path("output", "eac_cci_rf_pca.png"),
       plot = p_rf, width = 800 / 96, height = 600 / 96, dpi = 96, device = png)

month_data_compare <- month_data %>%
  select(trip_month, eac_cci) %>%
  left_join(month_data_rf, by = "trip_month") %>%
  pivot_longer(cols = c(eac_cci, eac_cci_rf),
               names_to = "index", values_to = "value")

p_compare <- ggplot(month_data_compare, aes(trip_month, value, colour = index)) +
  geom_line() +
  labs(x = "Month", y = "Index value",
       title = "EAC copepod composition index: Original vs RF + PCA")

ggsave(file.path("output", "eac_cci_comparison.png"),
       plot = p_compare, width = 800 / 96, height = 600 / 96, dpi = 96, device = png)

### --- Optional: Save other diagnostics from original method ---

samples_t <- samples %>%
  mutate(calc_seas_vars(sample_time)) %>%
  select(latitude, win_sum, aut_spr, sst)

png(file.path("output", "eac_cci_rda_envfit.png"), width = 680, height = 680, res = 96)
plot(rda_fit, display = "wa", choices = c(1, 2), scaling = "sites",
     main = format(terms(rda_fit)), cex.main = 1)
plot(envfit(rda_fit, samples_t))
invisible(dev.off())

png(file.path("output", "eac_cci_climatology_gam_terms.png"), width = 680, height = 680, res = 96)
plot(lm_fit, pages = 1)
invisible(dev.off())

p_clim <- ggplot(climatology, aes(sample_time, eac_cci)) +
  geom_line() +
  scale_x_datetime(date_breaks = "1 month", minor_breaks = NULL,
                   date_labels = "%b")
ggsave(file.path("output", "eac_cci_climatology.png"),
       plot = p_clim, width = 800 / 96, height = 600 / 96, dpi = 96, device = png)



month_data_rf_rda <- cci_data$month_data_rf_rda

month_data_t_rf_rda <- month_data_rf_rda %>%
  mutate(
    sample_time = trip_month + ddays(days_in_month(trip_month)) * 0.5,
    sample_time = floor_date(sample_time, unit = "day"),
    doy = yday(sample_time)
  ) %>%
  left_join(climatology %>% select(doy, eac_cci_clim = eac_cci), by = "doy") %>%
  mutate(anom_label = if_else(eac_cci_rf_rda >= eac_cci_clim, "Anomaly (+)", "Anomaly (-)"))

time_range <- range(month_data_t_rf_rda$sample_time)

clim_points_rf_rda <- tibble(sample_time = seq(time_range[1], time_range[2], by = "1 day")) %>%
  mutate(doy = pmin(yday(sample_time), 365)) %>%
  left_join(climatology %>% select(doy, eac_cci), by = "doy")

p_rf_rda <-
  ggplot(mapping = aes(sample_time, eac_cci_rf_rda)) +
  geom_segment(data = month_data_t_rf_rda,
               mapping = aes(xend = sample_time, yend = eac_cci_clim,
                             colour = anom_label, linetype = anom_label)) +
  geom_line(data = clim_points_rf_rda,
            mapping = aes(y = eac_cci, linetype = "Climatology", colour = "Climatology")) +
  geom_point(data = month_data_t_rf_rda, size = 2, shape = 21, fill = "gray") +
  scale_colour_manual(breaks = c("Anomaly (+)", "Anomaly (-)", "Climatology"),
                      values = c("red", "blue", "black")) +
  scale_linetype_manual(breaks = c("Anomaly (+)", "Anomaly (-)", "Climatology"),
                        values = c("solid", "dashed", "dotted")) +
  scale_x_datetime(date_breaks = "1 year", minor_breaks = NULL,
                   date_labels = "%Y") +
  labs(x = "Time",
       y = "EAC copepod composition index (RDA + RF)",
       colour = NULL,
       linetype = NULL,
       title = "EAC copepod composition index (RDA + RF, monthly average)")

ggsave(file.path("output", "eac_cci_rf_rda.png"),
       plot = p_rf_rda, width = 800 / 96, height = 600 / 96, dpi = 96, device = png)


climatology_rf_rda <- cci_data$climatology_rf_rda

p_clim_rf_rda <- ggplot(climatology_rf_rda, aes(doy, eac_cci_rf_rda)) +
  geom_line() +
  labs(x = "Day of year",
       y = "EAC CCI (RDA + RF)",
       title = "Climatology of EAC copepod composition index (RDA + RF)")

ggsave(file.path("output", "eac_cci_climatology_rf_rda.png"),
       plot = p_clim_rf_rda, width = 800 / 96, height = 600 / 96, dpi = 96)



# compare the three indexes

# Join and reshape all monthly data
comparison_data <- month_data %>%
  select(trip_month, eac_cci) %>%
  left_join(month_data_rf, by = "trip_month") %>%
  left_join(month_data_rf_rda, by = "trip_month") %>%
  pivot_longer(cols = c(eac_cci, eac_cci_rf, eac_cci_rf_rda),
               names_to = "index", values_to = "value")

p_compare_all <- ggplot(comparison_data, aes(trip_month, value, colour = index)) +
  geom_line(linewidth = 0.7) +
  labs(x = "Time",
       y = "Index value",
       title = "Comparison of EAC copepod composition indexes",
       colour = "Index type") +
  theme_minimal()

ggsave(file.path("output", "eac_cci_comparison_all.png"),
       plot = p_compare_all, width = 800 / 96, height = 600 / 96, dpi = 96)


# plotting comparisons and such

# Merge climatology to each monthly dataset
month_data_long <- bind_rows(
  month_data %>% mutate(index = "eac_cci"),
  month_data_rf %>% rename(eac_cci = eac_cci_rf) %>% mutate(index = "eac_cci_rf"),
  month_data_rf_rda %>% rename(eac_cci = eac_cci_rf_rda) %>% mutate(index = "eac_cci_rf_rda")
) %>%
  mutate(doy = yday(trip_month))

climatology_long <- bind_rows(
  climatology %>% mutate(index = "eac_cci", clim = eac_cci),
  climatology_rf %>% mutate(index = "eac_cci_rf", clim = eac_cci_rf),
  climatology_rf_rda %>% mutate(index = "eac_cci_rf_rda", clim = eac_cci_rf_rda)
)

# Join climatology
month_data_long <- month_data_long %>%
  left_join(climatology_long %>% select(doy, index, clim), by = c("doy", "index")) %>%
  mutate(anomaly_type = case_when(
    eac_cci > clim ~ "Positive anomaly",
    eac_cci < clim ~ "Negative anomaly",
    TRUE ~ NA_character_
  ))

p_anom <- ggplot(month_data_long, aes(trip_month, eac_cci)) +
  geom_segment(aes(xend = trip_month, yend = clim, colour = anomaly_type), linewidth = 0.5) +
  geom_line(aes(y = clim), linetype = "dotted", colour = "black") +
  geom_point(size = 1.2, shape = 21, fill = "white") +
  scale_colour_manual(values = c("Positive anomaly" = "red", "Negative anomaly" = "blue")) +
  labs(x = "Time", y = "Index value",
       title = "EAC copepod composition indexes with anomalies vs. climatology",
       colour = "Anomaly type") +
  facet_wrap(~ index, ncol = 1, scales = "free_y") +
  theme_minimal()

ggsave("output/eac_cci_all_anomalies.png", plot = p_anom,
       width = 900 / 96, height = 1200 / 96, dpi = 96)

p_facet <- ggplot(month_data_long, aes(trip_month, eac_cci)) +
  geom_line(aes(colour = index)) +
  labs(x = "Time", y = "Index value",
       title = "EAC copepod composition indexes (facet view)") +
  facet_wrap(~ index, ncol = 1, scales = "free_y") +
  theme_minimal()

ggsave("output/eac_cci_facet_plot.png", plot = p_facet,
       width = 900 / 96, height = 1200 / 96, dpi = 96)


# smooth trends

p_smooth <- ggplot(month_data_long, aes(trip_month, eac_cci, colour = index)) +
  geom_line(alpha = 0.4) +
  geom_smooth(method = "loess", span = 0.3, se = FALSE, linewidth = 1) +
  labs(x = "Time", y = "Index value",
       title = "Smoothed trends in EAC copepod composition indexes") +
  theme_minimal()

ggsave("output/eac_cci_smoothed_trends.png", plot = p_smooth,
       width = 900 / 96, height = 600 / 96, dpi = 96)
