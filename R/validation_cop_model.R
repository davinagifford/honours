### validation_cop_model.R
###
### compare output of index with data from the Copernicus Marine Service
###
### Created: 2025-08-19
### Author: Davina Gifford
### Last updated: 2025-09-09
### Edited by: Davina Gifford


library(tidyverse)
library(ncdf4)
library(tidync)
library(mgcv)
library(ggplot2)



# Use Copernicus data for validation --------------------------------------
# Global Ocean Physics Reanalysis Products


# Load the netCDF file

datafile <- "D:/HONOURS/DavinaG_2025_Honours/data/cmems_mod_glo_phy_my_0.083deg_P1M-m_1757126612458.nc"
if (!tolower(Sys.info()[["sysname"]]) == "sunos") {
  tnc <- tidync(datafile)
  print(tnc)
}

# explore variables
variables <- hyper_vars(tnc)
print(variables)

# extract velocity data
model_vel <- tnc %>%
  hyper_tibble(select_var = "vo")
print(model_vel)

raw_model_vel <- tnc %>%
  hyper_tibble(select_var = "vo")

raw_model_vel <- raw_model_vel %>% 
  filter(depth < 60) %>%
  filter(longitude > 153.9) %>% 
  filter(longitude < 155.0) 

# load the interim data

datafile2 <- "D:/HONOURS/DavinaG_2025_Honours/data/cmems_mod_glo_phy_myint_0.083deg_P1M-m_1757305199262.nc"
if (!tolower(Sys.info()[["sysname"]]) == "sunos") {
  tnc2 <- tidync(datafile2)
  print(tnc2)
}

# explore variables
variables <- hyper_vars(tnc)
print(variables)

# extract velocity data
model_vel2 <- tnc2 %>%
  hyper_tibble(select_var = "vo")
print(model_vel2)

raw_model_vel2 <- tnc2 %>%
  hyper_tibble(select_var = "vo")

raw_model_vel2 <- raw_model_vel2 %>% 
  filter(depth < 60) %>%
  filter(longitude > 153.9) %>% 
  filter(longitude < 155.0) 

model_vel$depth <- as.numeric(model_vel$depth)
model_vel2$depth <- as.numeric(model_vel2$depth)


model_vel2_50 <- model_vel2 %>%
  filter(depth < 60) %>%
  filter(longitude > 153.9) %>% 
  filter(longitude < 155.0)
  mutate(
    date = as.Date(sub("T.*", "", time)),
    year = year(date),
    month = month(date)
  ) %>%
  group_by(year, month) %>%
  summarise(mean_vel = mean(vo, na.rm = TRUE), .groups = "drop") %>% 
  mutate(date = as.Date(paste(year, month, "01", sep = "-"))) 


# Combine the two datasets

model_vel_full_50 <- bind_rows(model_vel_50, model_vel2_50)

# convert velocity to positive strength value
model_vel_full_50 <- model_vel_full_50 %>%
  mutate(mean_vel = abs(mean_vel))


# plot model velocity
p <- ggplot(model_vel_full_50, aes(x = date, y = mean_vel)) +
  geom_line(colour = "blue", linewidth = 1) +
  geom_point(size = 3, shape = 21, fill = "blue") +
  geom_smooth(method = "loess", se = TRUE, color = "blue", linetype = "dashed") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    x = "Time",
    y = expression("Monthly Mean model velocity (m/s"^{-1}*")"),
    title = "Copernicus Model Monthly Mean Velocity (50m depth)"
  )

ggsave(file.path("output", "model_velocity_50.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)



# 50m depth ---------------------------------------------------------------


# not the full copernicus data
model_vel$depth <- as.numeric(model_vel$depth)

model_vel_50 <- model_vel %>%
  filter(depth < 60) %>%
  filter(longitude > 153.9) %>% 
  filter(longitude < 155.0) %>% 
  mutate(
    date = as.Date(sub("T.*", "", time)),
    year = year(date),
    month = month(date)
  ) %>%
  group_by(year, month) %>%
  summarise(mean_vel = mean(vo, na.rm = TRUE), .groups = "drop") %>% 
  mutate(date = as.Date(paste(year, month, "01", sep = "-"))) 


# convert velocity to positive strength value
model_vel_50 <- model_vel_50 %>%
  mutate(mean_vel = abs(mean_vel))


# plot model velocity
p <- ggplot(model_vel_50, aes(x = date, y = mean_vel)) +
  geom_line(colour = "blue", linewidth = 1) +
  geom_point(size = 3, shape = 21, fill = "blue") +
  geom_smooth(method = "loess", se = TRUE, color = "blue", linetype = "dashed") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    x = "Time",
    y = expression("Monthly Mean model velocity (m/s"^{-1}*")"),
    title = "Copernicus Model Monthly Mean Velocity (50m depth)"
  )

ggsave(file.path("output", "model_velocity_50.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)



# Strength against climatology --------------------------------------------

# climatology for the 2011 - 2021 data

cop_model_data_50 <- model_vel_50 %>%
  mutate(
    doy = yday(date)
  )

lm_fit_strength_50 <- gam(mean_vel ~ s(doy, bs = "cc", k = 5),
                       data = cop_model_data_50,
                       method = "REML",
                       knots = list(doy = c(0, 365)))

mean_time_50 <- median(cop_model_data_50$date)
time0_50 <- min(cop_model_data_50$date)

climatology_str_50 <-
  tibble(date = floor_date(mean_time_50, unit = "year") + days(0:364),
         doy = yday(date),
         time_x = time_length(interval(time0_50, mean_time_50), unit = "day"))

clim_terms_str_50 <- predict(lm_fit_strength_50, type = "terms", newdata = climatology_str_50)

head(clim_terms_str_50)

climatology_str_50 <-
  climatology_str_50 %>%
  mutate(doy_eff = clim_terms_str_50[, "s(doy)"],
         intercept = coef(lm_fit_strength_50)["(Intercept)"],
         str_clim = doy_eff + intercept) %>%
  select(!c(doy_eff, intercept))

climatology_str_50 <- climatology_str_50 %>%
  mutate(date = as.Date(date))


month_clim_str_50 <- climatology_str_50 %>%
  mutate(month = month(date)) %>%
  group_by(month) %>%
  summarise(month_clim_str = mean(str_clim, na.rm = TRUE), .groups = "drop")

clim_compare_short_50 <- month_clim_str_50 %>%
  select(month, month_clim_str) %>%
  left_join(
    month_climatology,
    by = "month"
  )



# Pivot longer for easier plotting
clim_compare_long_50 <- clim_compare_short_50 %>%
  pivot_longer(cols = c(month_clim_str, month_clim),
               names_to = "climatology_type",
               values_to = "value")

# Calculate correlation


# Run Pearson correlation
correlation_50 <- cor(clim_compare_short_50$month_clim, clim_compare_short_50$month_clim_str, method = "pearson", use = "complete.obs")
cor_test_50 <- cor.test(clim_compare_short_50$month_clim, clim_compare_short_50$month_clim_str, method = "pearson", use = "complete.obs")

# Output results
print(paste("Pearson correlation:", round(correlation_50, 3)))
print(cor_test_50)


p_val_50 <- cor_test_50$p.value

# Run kendall correlation
correlation_50_k <- cor(clim_compare_short_50$month_clim, clim_compare_short_50$month_clim_str, method = "kendall", use = "complete.obs")
cor_test_50_k <- cor.test(clim_compare_short_50$month_clim, clim_compare_short_50$month_clim_str, method = "kendall", use = "complete.obs")

# Output results
print(paste("Kendall correlation:", round(correlation_50_k, 3)))
print(cor_test_50_k)


p_val_50_k <- cor_test_50_k$p.value

# Plot comparison
p <- ggplot(clim_compare_long_50, aes(x = month, y = value, color = climatology_type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3, shape = 21, fill = "white") +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_color_manual(values = c("month_clim_str" = "blue", "month_clim" = "red"),
                     labels = c("EAC CCI Climatology", "Current Strength Climatology")) +
  labs(
    x = "Month",
    y = "Climatology Value",
    color = "Climatology Type",
    title = "Comparison of Current Strength and EAC CCI Climatologies",
    subtitle = paste("Pearson correlation:", round(correlation_50, 3), "| p-value:", signif(p_val_50, 3))
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )


ggsave(file.path("output", "climatology-comparison_50m.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)


# plot a scatterplot of the two

ggplot(data = clim_compare_short_50) +
  geom_point(mapping = aes(x = month_clim, y = month_clim_str)) +
  geom_smooth(mapping = aes(x = month_clim, y = month_clim_str), method = "lm", se = TRUE, color = "blue", linetype = "dashed") +
  labs(
    x = "EAC CCI Climatology",
    y = "Current Strength Climatology",
    title = "Scatterplot of Current Strength vs EAC CCI Climatologies",
    subtitle = paste("Pearson correlation:", round(correlation_50, 3), "| p-value:", signif(p_val_50, 3))
  ) 





# climatology for full copernicus data-----------------

# calculate a climatology


cop_model_data_full_50 <- model_vel_full_50 %>%
  mutate(
    doy = yday(date)
  )

lm_fit_strength_full_50 <- gam(mean_vel ~ s(doy, bs = "cc", k = 5),
                        data = cop_model_data_full_50,
                        method = "REML",
                        knots = list(doy = c(0, 365)))

mean_time_full_50 <- median(cop_model_data_full_50$date)
time0_full_50 <- min(cop_model_data_full_50$date)

climatology_str_full_50 <-
  tibble(date = floor_date(mean_time_full_50, unit = "year") + days(0:364),
         doy = yday(date),
         time_x = time_length(interval(time0_full_50, mean_time_full_50), unit = "day"))

clim_terms_str_full_50 <- predict(lm_fit_strength_full_50, type = "terms", newdata = climatology_str_full_50)

head(clim_terms_str_full_50)

climatology_str_full_50 <-
  climatology_str_full_50 %>%
  mutate(doy_eff = clim_terms_str_full_50[, "s(doy)"],
         intercept = coef(lm_fit_strength_full_50)["(Intercept)"],
         str_clim = doy_eff + intercept) %>%
  select(!c(doy_eff, intercept))

climatology_str_full_50 <- climatology_str_full_50 %>%
  mutate(date = as.Date(date))


month_clim_str_full_50 <- climatology_str_full_50 %>%
  mutate(month = month(date)) %>%
  group_by(month) %>%
  summarise(month_clim_str = mean(str_clim, na.rm = TRUE), .groups = "drop")

clim_compare_short_full_50 <- month_clim_str_full_50 %>%
  select(month, month_clim_str) %>%
  left_join(
    month_climatology,
    by = "month"
  )



# Pivot longer for easier plotting
clim_compare_long_full_50 <- clim_compare_short_full_50 %>%
  pivot_longer(cols = c(month_clim_str, month_clim),
               names_to = "climatology_type",
               values_to = "value")

# Calculate correlation


# Run Pearson correlation
correlation_full_50 <- cor(clim_compare_short_full_50$month_clim, clim_compare_short_full_50$month_clim_str, method = "pearson", use = "complete.obs")
cor_test_full_50 <- cor.test(clim_compare_short_full_50$month_clim, clim_compare_short_full_50$month_clim_str, method = "pearson", use = "complete.obs")

# Output results
print(paste("Pearson correlation:", round(correlation_full_50, 3)))
print(cor_test_full_50)


p_val_full_50 <- cor_test_full_50$p.value

# Run kendall correlation
correlation_full_50_k <- cor(clim_compare_short_full_50$month_clim, clim_compare_short_full_50$month_clim_str, method = "kendall", use = "complete.obs")
cor_test_full_50_k <- cor.test(clim_compare_short_full_50$month_clim, clim_compare_short_full_50$month_clim_str, method = "kendall", use = "complete.obs")

# Output results
print(paste("Kendall correlation:", round(correlation_full_50_k, 3)))
print(cor_test_full_50_k)


p_val_full_50_k <- cor_test_full_50_k$p.value

# Plot comparison
p <- ggplot(clim_compare_long_full_50, aes(x = month, y = value, color = climatology_type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3, shape = 21, fill = "white") +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_color_manual(values = c("month_clim_str" = "blue", "month_clim" = "red"),
                     labels = c("EAC CCI Climatology", "Current Strength Climatology")) +
  labs(
    x = "Month",
    y = "Climatology Value",
    color = "Climatology Type",
    title = "Comparison of Current Strength and EAC CCI Climatologies - (50 m depth)",
    subtitle = paste("Pearson correlation:", round(correlation_full_50, 3), "| p-value:", signif(p_val_full_50, 3))
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )


ggsave(file.path("output", "climatology-comparison_full_50.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)


# plot a scatterplot of the two

ggplot(data = clim_compare_short_full_50) +
  geom_point(mapping = aes(x = month_clim, y = month_clim_str)) +
  geom_smooth(mapping = aes(x = month_clim, y = month_clim_str), method = "lm", se = TRUE, color = "blue", linetype = "dashed") +
  labs(
    x = "EAC CCI Climatology",
    y = "Current Strength Climatology",
    title = "Scatterplot of Current Strength vs EAC CCI Climatologies",
    subtitle = paste("Pearson correlation:", round(correlation_full_50, 3), "| p-value:", signif(p_val_full_50, 3))
  ) 




# Residuals  --------------------------------------------------------------
# linear regression of the two climatologies with the copernicus data from 2011 to 2021

# linear regression of the two climatologies

# Fit a linear model
cop_model_50 <- lm(month_clim ~ month_clim_str, data = clim_compare_short_50)

# Summary of the model
summary(cop_model_50)



residuals_cop_50 <- resid(cop_model_50)
predicted_cop_50 <- fitted(cop_model_50)



plot(clim_compare_short_50$month_clim_str, cop_model_50$residuals,
     main = "Residuals vs Fitted Values",
     xlab = "Fitted Values",
     ylab = "Residuals",
     pch = 19, col = "darkgreen")
abline(h = 0, col = "red", lwd = 2)





# linear regression of the two climatologies with full copernicus data

# Fit a linear model
cop_model_full_50 <- lm(month_clim ~ month_clim_str, data = clim_compare_short_full_50)

# Summary of the model
summary(cop_model_full_50)



residuals_cop_full_50 <- resid(cop_model_full_50)
predicted_cop_full_50 <- fitted(cop_model_full_50)



plot(clim_compare_short_full_50$month_clim_str, cop_model_full_50$residuals,
     main = "Residuals vs Fitted Values",
     xlab = "Fitted Values",
     ylab = "Residuals",
     pch = 19, col = "darkgreen")
abline(h = 0, col = "red", lwd = 2)





# Strength against climatology --------------------------------------------

time_range_full_50 <- cop_model_data_full_50 %>%
  reframe(time_range = range(date)) %>%
  pull()


clim_points_full_50 <- tibble(date = seq(time_range_full_50[1], time_range_full_50[2], by = "1 day")) %>%
  mutate(doy = pmin(yday(date), 365)) %>%
  left_join(climatology_str_full_50 %>% select(doy, str_clim), by = "doy")

clim_points_full_50 <- clim_points_full_50 %>% 
  mutate(date = as.POSIXct(date))





strength_with_clim_full_50 <- cop_model_data_full_50 %>%
  mutate(doy = yday(date)) %>%
  left_join(climatology_str_full_50 %>% select(doy, str_clim), by = "doy") %>%
  mutate(
    vel_anomaly = mean_vel - str_clim,
    anom_label = case_when(
      vel_anomaly > 0 ~ "Anomaly (+)",
      TRUE ~ "Anomaly (-)")
  )


strength_with_clim_full_50 <- strength_with_clim_full_50 %>%
  mutate(date = as.POSIXct(date))


str_anom_full_50 <- strength_with_clim_full_50 %>% 
  select(vel_anomaly)

# Plot strength with climatology
p <- ggplot() +
  # Segments for anomalies
  geom_segment(
    data = strength_with_clim_full_50,
    aes(x = date, xend = date, y = str_clim, yend = mean_vel,
        colour = anom_label, linetype = anom_label),
    linewidth = 1.2
  ) +
  # Climatology line
  geom_line(
    data = clim_points_full_50,
    aes(x = date, y = str_clim, linetype = "Climatology", colour = "Climatology"),
    linewidth = 1.1
  ) +
  # Points for observed mean_vcur
  geom_point(
    data = strength_with_clim_full_50,
    aes(x = date, y = mean_vel, fill = anom_label, colour = anom_label),
    size = 2.5, shape = 21
  ) +
  scale_colour_manual(
    breaks = c("Anomaly (+)", "Anomaly (-)", "Climatology"),
    values = c("Anomaly (+)" = "red", "Anomaly (-)" = "blue", "Climatology" = "black")
  ) +
  scale_fill_manual(
    breaks = c("Anomaly (+)", "Anomaly (-)"),
    values = c("Anomaly (+)" = "red", "Anomaly (-)" = "blue"),
    guide = "none"
  ) +
  scale_linetype_manual(
    breaks = c("Anomaly (+)", "Anomaly (-)", "Climatology"),
    values = c("solid", "solid", "solid")
  ) +
  scale_x_datetime(date_breaks = "1 year", minor_breaks = NULL, date_labels = "%Y") +
  #coord_cartesian(ylim = c(-0.3, 0.4)) +
  labs(
    x = "Time",
    y = expression("Monthly Mean strength (m/s"^{-1}*")"),
    colour = NULL,
    fill = NULL,
    linetype = NULL,
    title = "Current Strength Climatology (monthly average) (50m depth)"
  )

ggsave(file.path("output", "strength-with-clim_full_50.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)


# data wrangling

# data wrangling

raw_strength_data <- model_vel_full_50 %>% 
  select(date, mean_vel)

raw_index_data <- month_data_t %>% 
  select(trip_month, eac_cci)

raw_index_data <- raw_index_data %>% 
  rename(date = trip_month)

raw_index_data <- raw_index_data %>% 
  mutate(
    date = as.Date(sub("T.*", "", trip_month)),
    year = year(date),
    month = month(date)
  ) 
