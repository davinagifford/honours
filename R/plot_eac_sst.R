### plot_eac_sst.R
###
### Plots EAC monthly nominal SST.
###
### Created: 2023-07-21
### Author: Wayne A. Rochester
### Last updated: 2025-04-28
### Edited by: Davina Gifford

library(lubridate)
library(mgcv)
library(dplyr)
library(tidyr)
library(ggplot2)

sst_data <- readRDS(file.path("var", "eac_sst.rds"))
cci_data <- readRDS(file.path("var", "eac_cci.rds"))

samples_new <- sst_data$samples
lm_fit <- sst_data$lm_fit
climatology <- sst_data$climatology
month_data <- sst_data$month_data

time_range<-
    cci_data$month_data %>%
    mutate(sample_time = trip_month + ddays(days_in_month(trip_month)) * 0.5,
           sample_time = floor_date(sample_time, unit = "day")) %>%
    reframe(time_range = range(sample_time)) %>%
    pull()
time_range <- time_range + days(c(-90, 90))

month_data_t <-
    month_data %>%
    mutate(sample_time = sample_month +
               ddays(days_in_month(sample_month)) * 0.5,
           sample_time = floor_date(sample_time, unit = "day"),
           doy = yday(sample_time)) %>%
    left_join(climatology %>% select(doy, eac_sst_clim = eac_sst),
              by = "doy") %>%
    mutate(anom_label = if_else(eac_sst >= eac_sst_clim,
                                "Anomaly (+)",
                                "Anomaly (-)")) %>%
    filter(sample_time >= time_range[1], sample_time <= time_range[2])

time_range <-
    month_data_t %>%
    reframe(time_range = range(sample_time)) %>%
    pull()
xlim <- time_range + diff(range(time_range)) * c(-0.03, 0.03)
ylim <- month_data_t %>% reframe(ylim = range(eac_sst)) %>% pull()
ylim <- ylim + c(-1, 1) * diff(ylim) * 0.05

clim_points <-
    tibble(sample_time = seq(time_range[1], time_range[2], by = "1 day")) %>%
    mutate(doy = pmin(yday(sample_time), 365)) %>%
    left_join(climatology %>% select(doy, eac_sst), by = "doy")

p <-
    ggplot(mapping = aes(sample_time, eac_sst)) +
    geom_segment(data = month_data_t,
                 mapping = aes(xend = sample_time, yend = eac_sst_clim,
                               colour = anom_label, linetype = anom_label)) +
    geom_line(data = clim_points,
              mapping = aes(linetype = "Climatology",
                            colour = "Climatology")) +
    geom_point(data = month_data_t, size = 2, shape = 21, fill = "gray") +
    scale_colour_manual(breaks = c("Anomaly (+)", "Anomaly (-)",
                                   "Climatology"),
                        values = c("red", "blue", "black")) +
    scale_linetype_manual(breaks = c("Anomaly (+)", "Anomaly (-)",
                                     "Climatology"),
                          values = c("solid", "dashed", "dotted")) +
    scale_x_datetime(date_breaks = "1 year", minor_breaks = NULL,
                     date_labels = "%Y") +
    coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(x = "Time",
         y = expression(paste("SST ("*degree, "C)")),
         colour = NULL,
         linetype = NULL,
         title = "EAC nominal SST (monthly average)")
ggsave(file.path("output", "eac_sst.png"),
       plot = p, width = 800 / 96, height = 600 / 96, dpi = 96,
       device = png)


png(file.path("output", "eac_sst_climatology_gam_terms.png"),
    width = 680, height = 680, res = 96)
plot(lm_fit, pages = 1)
invisible(dev.off())


p <-
    ggplot(climatology, aes(sample_time, eac_sst)) +
    geom_line() +
    scale_x_datetime(date_breaks = "1 month", minor_breaks = NULL,
                     date_labels = "%b")
ggsave(file.path("output", "eac_sst_climatology.png"),
       plot = p, width = 800 / 96, height = 600 / 96, dpi = 96,
       device = png)


moy_data <-
    month_data %>%
    mutate(moy = month(sample_month)) %>%
    group_by(moy) %>%
    summarise(eac_sst = mean(eac_sst), n = n(), .groups = "drop")

clim_moy_data <-
    climatology %>%
    mutate(moy = month(sample_time)) %>%
    group_by(moy) %>%
    summarise(eac_sst = mean(eac_sst), n = n(), .groups = "drop")

clim_check_data <-
    month_data %>%
    mutate(moy = month(sample_month), src_code = "obs") %>%
    select(src_code, moy, eac_sst) %>%
    bind_rows(moy_data %>%
              mutate(src_code = "avg") %>%
              select(src_code, moy, eac_sst)) %>%
    bind_rows(clim_moy_data %>%
              mutate(src_code = "clim") %>%
              select(src_code, moy, eac_sst))

p <-
    ggplot(clim_check_data, aes(moy, eac_sst)) +
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
ggsave(file.path("output", "eac_sst_climatology_check.png"),
       plot = p, width = 800 / 96, height = 600 / 96, dpi = 96,
       device = png)
