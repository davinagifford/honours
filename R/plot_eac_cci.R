### plot_eac_cci.R
###
### Plots the EAC copepod composition index.
###
### Created: 2023-07-21
### Author: Wayne A. Rochester
### Last updated: 2026-01-29
### Edited by: Davina Gifford

library(lubridate)
library(vegan)
library(mgcv)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(scales)     # for colour scales
library(ggrepel)    # for readable labels
library(ggvegan)

install.packages("remotes")
remotes::install_github("gavinsimpson/ggvegan")


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


month_data_t <-
    month_data %>%
    mutate(sample_time = trip_month + ddays(days_in_month(trip_month)) * 0.5,
           sample_time = floor_date(sample_time, unit = "day"),
           doy = yday(sample_time)) %>%
    left_join(climatology %>% select(doy, eac_cci_clim = eac_cci),
              by = "doy") %>%
    mutate(anom_label = if_else(eac_cci >= eac_cci_clim,
                                "Anomaly (+)",
                                "Anomaly (-)"))

time_range <-
    month_data_t %>%
    reframe(time_range = range(sample_time)) %>%
    pull()
clim_points <-
    tibble(sample_time = seq(time_range[1], time_range[2], by = "1 day")) %>%
    mutate(doy = pmin(yday(sample_time), 365)) %>%
    left_join(climatology %>% select(doy, eac_cci), by = "doy")

# Plot index without climatology
p <-
  ggplot(month_data_t, aes(x = sample_time, y = eac_cci)) +
  geom_line(linewidth = 1, colour = "black") +
  geom_point(size = 3, shape = 21, fill = "black") +
  scale_x_datetime(date_breaks = "1 year", minor_breaks = NULL,
                   date_labels = "%Y") +
  geom_smooth(method = "loess", se = TRUE, color = "black", linetype = "dashed") +
  theme(axis.title = element_text(size = 14), 
        plot.title = element_text(size = 18, face = "bold")) +
  labs(x = "Time",
       y = "EAC copepod composition index",
       title = "EAC copepod composition index (monthly average)") 

ggsave(file.path("output", "eac_cci_noclim.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)

# Plot index with climatology
p <-
    ggplot(mapping = aes(sample_time, eac_cci)) +
    geom_segment(data = month_data_t,
                 mapping = aes(xend = sample_time, yend = eac_cci_clim,
                               colour = anom_label, linetype = anom_label),
                 linewidth = 1.5) +
    geom_line(data = clim_points,
              mapping = aes(linetype = "Climatology",
                            colour = "Climatology"),
              linewidth = 1) +
    geom_segment(data = month_data_t,
               mapping = aes(xend = sample_time, yend = eac_cci_clim,
                             colour = anom_label, linetype = anom_label),
               linewidth = 1.5) +
    geom_point(data = month_data_t, size = 3, shape = 21, fill = "black") +
    geom_point(data = month_data_t, size = 2.5, shape = 21, fill = "black") +
    scale_colour_manual(breaks = c("Anomaly (+)", "Anomaly (-)",
                                   "Climatology"),
                        values = c("red", "blue", "black")) +
    scale_linetype_manual(breaks = c("Anomaly (+)", "Anomaly (-)",
                                     "Climatology"),
                          values = c("solid", "solid", "solid")) +
    scale_x_datetime(date_breaks = "1 year", minor_breaks = NULL,
                     date_labels = "%Y") +
    theme(axis.title = element_text(size = 14), 
          plot.title = element_text(size = 18, face = "bold")) +
    labs(x = "Time",
         y = "EAC copepod composition index",
         colour = NULL,
         linetype = NULL,
         title = "EAC copepod composition index (monthly average)",
         )
ggsave(file.path("output", "eac_cci.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)


samples_t <-
    samples %>%
    mutate(calc_seas_vars(sample_time)) %>%
    select(latitude, win_sum, aut_spr, sst)

png(file.path("output", "eac_cci_rda_envfit_new.png"),
    width = 680, height = 680, res = 96)
plot(rda_fit, display = "wa", choices = c(1, 2), scaling = "sites",
     main = format(terms(rda_fit)), cex.main = 1)
plot(envfit(rda_fit, samples_t))

invisible(dev.off())


png(file.path("output", "eac_cci_climatology_gam_terms_new.png"),
    width = 680, height = 680, res = 96)
plot(lm_fit, pages = 1)
invisible(dev.off())


p <-
    ggplot(climatology, aes(sample_time, eac_cci)) +
    geom_line() +
    scale_x_datetime(date_breaks = "1 month", minor_breaks = NULL,
                     date_labels = "%b")
ggsave(file.path("output", "eac_cci_climatology_new.png"),
       plot = p, width = 800 / 96, height = 600 / 96, dpi = 96,
       device = png)

# check the climatology to ensure it is logical
moy_data <-
    month_data %>%
    mutate(moy = month(trip_month)) %>%
    group_by(moy) %>%
    summarise(eac_cci = mean(eac_cci), n = n(), .groups = "drop")

clim_moy_data <-
    climatology %>%
    mutate(moy = month(sample_time)) %>%
    group_by(moy) %>%
    summarise(eac_cci = mean(eac_cci), n = n(), .groups = "drop")

clim_check_data <-
    month_data %>%
    mutate(moy = month(trip_month), src_code = "obs") %>%
    select(src_code, moy, eac_cci) %>%
    bind_rows(moy_data %>%
              mutate(src_code = "avg") %>%
              select(src_code, moy, eac_cci)) %>%
    bind_rows(clim_moy_data %>%
              mutate(src_code = "clim") %>%
              select(src_code, moy, eac_cci))


# check the climatology
p <-
    ggplot(clim_check_data, aes(moy, eac_cci)) +
    geom_line(data = clim_check_data %>% filter(src_code == "clim")) +
    geom_point(mapping = aes(shape = src_code, size = src_code,
                             colour = src_code, fill = src_code)) +
    scale_x_continuous(breaks = scales::breaks_width(1),
                       minor_breaks = NULL) +
    scale_shape_manual(breaks = c("obs", "avg", "clim"),
                       values = c(21, 1, 19)) +
    scale_size_manual(breaks = c("obs", "avg", "clim"),
                      values = c(3, 5, 2)) +
    scale_colour_manual(breaks = c("obs", "avg", "clim"),
                        values = c("black", "red", "black")) +
    scale_fill_manual(breaks = c("obs", "avg", "clim"),
                      values = c("gray", "red", "black"))
ggsave(file.path("output", "eac_cci_climatology_check_new.png"),
       plot = p, width = 800 / 96, height = 600 / 96, dpi = 96,
       device = png)


# graphs for thesis

plot_predictions(cci_str_model3, condition = "mean_strength", points= 0.25) + 
  labs(
    x = "mean_vel"
  )

plot_predictions(cci_str_model3, condition = "month_clim_str", points= 0.25) +
  labs(
    x = "vel_clim"
  )


plot_predictions(cci_str_model_n3, condition = "mean_strength", points= 0.25) + 
  labs(
    x = "mean_vel"
  )

plot_predictions(cci_str_model_n3, condition = "month_clim_str", points= 0.25)+
  labs(
    x = "vel_clim"
  )


plot_predictions(cci_str_model_s3, condition = "mean_strength", points= 0.25) +
  labs(
    x = "mean_vel"
  )

plot_predictions(cci_str_model_s3, condition = "month_clim_str", points= 0.25) +
  labs(
    x = "vel_clim"
  )


