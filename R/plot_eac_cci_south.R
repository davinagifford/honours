### plot_eac_cci.R
###
### Plots the EAC copepod composition index for the region North of the EAC Separation Zone.
###
### Created: 2023-07-21
### Author: Wayne A. Rochester
### Last updated: 2025-07-02
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


cci_data <- readRDS(file.path("var", "eac_cci_south.rds"))

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

p <-
    ggplot(mapping = aes(sample_time, eac_cci)) +
    geom_segment(data = month_data_t,
                 mapping = aes(xend = sample_time, yend = eac_cci_clim,
                               colour = anom_label, linetype = anom_label),
                 linewidth = 1.2) +
    geom_line(data = clim_points,
              mapping = aes(linetype = "Climatology",
                            colour = "Climatology"),
              linewidth = 1.1) +
    geom_point(data = month_data_t, size = 2.5, shape = 21, fill = "black") +
    scale_colour_manual(breaks = c("Anomaly (+)", "Anomaly (-)",
                                   "Climatology"),
                        values = c("red", "blue", "black")) +
    scale_linetype_manual(breaks = c("Anomaly (+)", "Anomaly (-)",
                                     "Climatology"),
                          values = c("solid", "solid", "solid")) +
    scale_x_datetime(date_breaks = "1 year", minor_breaks = NULL,
                     date_labels = "%Y") +
    labs(x = "Time",
         y = "EAC copepod composition index",
         colour = NULL,
         linetype = NULL,
         title = "EAC copepod composition index south of Separation Zone (monthly average)")
ggsave(file.path("output", "eac_cci_south.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)


samples_t <-
    samples %>%
    mutate(calc_seas_vars(sample_time)) %>%
    select(latitude, win_sum, aut_spr, sst)

png(file.path("output", "eac_cci_rda_envfit_south.png"),
    width = 680, height = 680, res = 96)
plot(rda_fit, display = "wa", choices = c(1, 2), scaling = "sites",
     main = format(terms(rda_fit)), cex.main = 1)
plot(envfit(rda_fit, samples_t))

invisible(dev.off())


png(file.path("output", "eac_cci_climatology_gam_terms_south.png"),
    width = 680, height = 680, res = 96)
plot(lm_fit, pages = 1)
invisible(dev.off())


p <-
    ggplot(climatology, aes(sample_time, eac_cci)) +
    geom_line() +
    scale_x_datetime(date_breaks = "1 month", minor_breaks = NULL,
                     date_labels = "%b")
ggsave(file.path("output", "eac_cci_climatology_south.png"),
       plot = p, width = 800 / 96, height = 600 / 96, dpi = 96,
       device = png)


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
ggsave(file.path("output", "eac_cci_climatology_check_south.png"),
       plot = p, width = 800 / 96, height = 600 / 96, dpi = 96,
       device = png)



# make it pretty
library(vegan)
library(ggplot2)

# 1. Fit envfit
fit_env <- envfit(rda_fit, samples_t)

# 2. Extract site scores
site_scores <- scores(rda_fit, display = "wa", choices = c(1, 2), scaling = "sites")
site_data <- as.data.frame(site_scores)
site_data$sst <- samples_t$sst

# 3. Extract and scale vectors

# Extract and scale vectors properly
vecs <- scores(fit_env, display = "vectors", scaling = "sites")  # This returns a matrix

# Convert to data frame and scale by vector length (magnitude)
vecs_df <- as.data.frame(vecs)
vecs_df$Variable <- rownames(vecs_df)

# Optionally scale vectors for visibility
vecs_df$RDA1 <- vecs_df$RDA1 * 0.3
vecs_df$PC1 <- vecs_df$PC1 * 0.23


# 4. Plot with ggplot2
ggplot(site_data, aes(x = RDA1, y = PC1, color = sst)) +
  geom_point(size = 3) +
  scale_color_viridis_c() +
  geom_segment(data = vecs_df,
               aes(x = 0, y = 0, xend = RDA1, yend = PC1),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black") +
  geom_text(data = vecs_df,
            aes(x = RDA1 * 1.1, y = PC1 * 1.1, label = Variable),
            color = "black", size = 4) +
  labs(title = format(terms(rda_fit)),
       x = "RDA1", y = "PC1") +
  theme_minimal()
