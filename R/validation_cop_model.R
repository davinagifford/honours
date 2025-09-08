### validation.R
###
### compare output of index with the original index
###
### Created: 2025-08-19
### Author: Davina Gifford
### Last updated: 2025-08-19
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
  filter(depth < 1000) %>%
  filter(longitude > 153.9) %>% 
  filter(longitude < 155.0)



model_vel <- model_vel %>%
  filter(depth < 1000) %>%
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

# plot model velocity
p <- ggplot(model_vel, aes(x = date, y = mean_vel)) +
  geom_line(colour = "blue", linewidth = 1) +
  geom_point(size = 3, shape = 21, fill = "blue") +
  geom_smooth(method = "loess", se = TRUE, color = "blue", linetype = "dashed") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    x = "Time",
    y = expression("Monthly Mean model velocity (m/s"^{-1}*")"),
    title = "Copernicus Model Monthly Mean Velocity (100m depth)"
  )

ggsave(file.path("output", "model_velocity.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)

# convert velocity to positive strength value
model_vel <- model_vel %>%
  mutate(mean_vel = abs(mean_vel))


# compare model velocity with observed strength (data_tbl)

model_str <- model_vel %>%
  left_join(data_tbl %>% select(date, mean_vcur), by = "date")
model_str <- model_str %>%
  drop_na()
# calculate correlation

correlation_model_str <- cor(model_str$mean_vel, model_str$mean_vcur, method = "pearson", use = "complete.obs")
print(paste("Pearson correlation between model velocity and observed strength:", correlation_model_str))

# compare with the index

model_vel_index <- model_vel %>%
  left_join(anomaly %>% select(date = trip_month, eac_cci = observed_eac_cci), by = "date")

model_vel_index <- model_vel_index %>%
  mutate(date = as.Date(date))


strength_scale <- max(model_vel_index$eac_cci, na.rm = TRUE) / max(abs(model_vel_index$mean_vel), na.rm = TRUE)

model_vel_index <- model_vel_index %>%
  mutate(vel_scale = mean_vel * strength_scale)

p <- ggplot(model_vel_index, aes(x = date)) +
  geom_line(aes(y = eac_cci), colour = "black", linewidth = 1) +
  geom_point(aes(y = eac_cci), size = 2, shape = 21, fill = "black") +
  geom_smooth(aes(y = eac_cci), method = "loess", se = TRUE, color = "black", linetype = "dashed") +
  geom_line(aes(y = mean_vel), colour = "blue", linewidth = 1) +
  geom_point(aes(y = mean_vel), size = 2, shape = 21, fill = "blue") +
  geom_smooth(aes(y = mean_vel), method = "loess", se = TRUE, color = "blue", linetype = "dashed") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(name = "EAC copepod composition index",
                     sec.axis = sec_axis(~ .,  name = "Monthy mean strength (Modelled) (m/s)")) +
  theme(
    axis.title.y = element_text(size = 14, color = "black"),
    axis.title.y.right = element_text(size = 14, color = "blue"),
    plot.title = element_text(size = 18, face = "bold")
  ) +
  labs(
    x = "Time",
    title = "EAC copepod composition index & Monthly mean strength (Copernicus Model data) ")



ggsave(file.path("output", "eac_cci_model_str_combined.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)

# test similarity between cci and strength values

model_vel_index <- model_vel_index %>% 
  mutate(
    year = year(date),
    month = month(date)
  ) 

model_str_combined_data_forcorr <- model_vel_index %>%
  select(month, eac_cci, mean_vel) %>%
  drop_na()

correlation_str <- cor(model_str_combined_data_forcorr$eac_cci, model_str_combined_data_forcorr$mean_vel, method = "pearson")
print(paste("Pearson correlation between EAC CCI and strength:", correlation_str)) 

# anova

anova_result_str <- aov(eac_cci ~ mean_vel, data = model_str_combined_data_forcorr)

summary(anova_result_str)

# linear regression of EAC CCI against strength

# Fit a linear model
model_str2 <- lm(eac_cci ~ mean_vel, data = model_vel_index)

# Summary of the model
summary(model_str2)



# scatterplot of model velocity against index 

ggplot(model_vel_index, aes(x = mean_vel, y = eac_cci)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(title = "Linear Regression of EAC CCI on Modelled Strength",
       x = expression("Monthly Mean strength (m/s"^{-1}*")"),
       y = "EAC CCI ") +
  theme_minimal()



# Strength against climatology --------------------------------------------

time_range <- current_data %>%
  reframe(time_range = range(date)) %>%
  pull()


clim_points <- tibble(date = seq(time_range[1], time_range[2], by = "1 day")) %>%
  mutate(doy = pmin(yday(date), 365)) %>%
  left_join(climatology_str %>% select(doy, str_clim), by = "doy")

clim_points <- clim_points %>% 
  mutate(date = as.POSIXct(date))





strength_with_clim <- current_data %>%
  mutate(doy = yday(date)) %>%
  left_join(climatology_str %>% select(doy, str_clim), by = "doy") %>%
  mutate(
    vcur_anomaly = mean_vcur - str_clim,
    anom_label = case_when(
      vcur_anomaly > 0 ~ "Anomaly (+)",
      TRUE ~ "Anomaly (-)")
  )


strength_with_clim <- strength_with_clim %>%
  mutate(date = as.POSIXct(date))


str_anom <- strength_with_clim %>% 
  select(vcur_anomaly)

# Plot strength with climatology
p <- ggplot() +
  # Segments for anomalies
  geom_segment(
    data = strength_with_clim,
    aes(x = date, xend = date, y = str_clim, yend = mean_vcur,
        colour = anom_label, linetype = anom_label),
    linewidth = 1.2
  ) +
  # Climatology line
  geom_line(
    data = clim_points,
    aes(x = date, y = str_clim, linetype = "Climatology", colour = "Climatology"),
    linewidth = 1.1
  ) +
  # Points for observed mean_vcur
  geom_point(
    data = strength_with_clim,
    aes(x = date, y = mean_vcur, fill = anom_label, colour = anom_label),
    size = 2.5, shape = 21
  ) +
  scale_colour_manual(
    breaks = c("Anomaly (+)", "Anomaly (-)", "Climatology"),
    values = c("Anomaly (+)" = "red", "Anomaly (-)" = "blue", "Climatology" = "black")
  ) +
  scale_fill_manual(
    breaks = c("Anomaly (+)", "Anomaly (-)"),
    values = c("Anomaly (+)" = "red", "Anomaly (-)" = "blue")
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
    title = "Current Strength Climatology (monthly average)"
  )

ggsave(file.path("output", "strength-with-clim.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)

# climatology for copernicus data-----------------

# calculate a climatology


cop_model_data <- model_vel %>%
  mutate(
    doy = yday(date)
  )

lm_fit_strength2 <- gam(mean_vel ~ s(doy, bs = "cc", k = 5),
                       data = cop_model_data,
                       method = "REML",
                       knots = list(doy = c(0, 365)))

mean_time <- median(cop_model_data$date)
time0 <- min(cop_model_data$date)

climatology_str2 <-
  tibble(date = floor_date(mean_time, unit = "year") + days(0:364),
         doy = yday(date),
         time_x = time_length(interval(time0, mean_time), unit = "day"))

clim_terms_str2 <- predict(lm_fit_strength2, type = "terms", newdata = climatology_str2)

head(clim_terms_str2)

climatology_str2 <-
  climatology_str2 %>%
  mutate(doy_eff = clim_terms_str2[, "s(doy)"],
         intercept = coef(lm_fit_strength2)["(Intercept)"],
         str_clim = doy_eff + intercept) %>%
  select(!c(doy_eff, intercept))

climatology_str2 <- climatology_str2 %>%
  mutate(date = as.Date(date))


month_clim_str2 <- climatology_str2 %>%
  mutate(month = month(date)) %>%
  group_by(month) %>%
  summarise(month_clim_str = mean(str_clim, na.rm = TRUE), .groups = "drop")

clim_compare_short2 <- month_clim_str2 %>%
  select(month, month_clim_str) %>%
  left_join(
    month_climatology,
    by = "month"
  )



# Pivot longer for easier plotting
clim_compare_long2 <- clim_compare_short2 %>%
  pivot_longer(cols = c(month_clim_str, month_clim),
               names_to = "climatology_type",
               values_to = "value")

# Calculate correlation


# Run Pearson correlation
correlation2 <- cor(clim_compare_short2$month_clim, clim_compare_short2$month_clim_str, method = "pearson", use = "complete.obs")
cor_test2 <- cor.test(clim_compare_short2$month_clim, clim_compare_short2$month_clim_str, method = "pearson", use = "complete.obs")

# Output results
print(paste("Pearson correlation:", round(correlation2, 3)))
print(cor_test2)


p_val2 <- cor_test2$p.value



# Plot comparison
p <- ggplot(clim_compare_long2, aes(x = month, y = value, color = climatology_type)) +
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
    subtitle = paste("Pearson correlation:", round(correlation2, 3), "| p-value:", signif(p_val2, 3))
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )


ggsave(file.path("output", "climatology-comparison2.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)


# plot a scatterplot of the two

ggplot(data = clim_compare_short2) +
  geom_point(mapping = aes(x = month_clim, y = month_clim_str)) +
  geom_smooth(mapping = aes(x = month_clim, y = month_clim_str), method = "lm", se = TRUE, color = "blue", linetype = "dashed") +
  labs(
    x = "EAC CCI Climatology",
    y = "Current Strength Climatology",
    title = "Scatterplot of Current Strength vs EAC CCI Climatologies",
    subtitle = paste("Pearson correlation:", round(correlation2, 3), "| p-value:", signif(p_val2, 3))
  ) 


# correlate the residuals

# Step 1: Calculate residuals (deviations from monthly means)
clim_compare_short2 <- clim_compare_short2 %>%
  mutate(
    residual_eac_cci = month_clim - mean(month_clim, na.rm = TRUE),
    residual_strength = month_clim_str - mean(month_clim_str, na.rm = TRUE)
  )

# Step 2: Correlate the residuals
residual_correlation2 <- cor(clim_compare_short2$residual_eac_cci, clim_compare_short2$residual_strength, method = "pearson", use = "complete.obs")
residual_cor_test2 <- cor.test(clim_compare_short2$residual_eac_cci, clim_compare_short2$residual_strength)

# Step 3: Print results
print(paste("Residual correlation:", round(residual_correlation2, 3)))
print(residual_cor_test2)

p_val3 <- cor_test2$p.value

# Scatterplot of residuals
ggplot(clim_compare_short2, aes(x = residual_eac_cci, y = residual_strength)) +
  geom_point(color = "blue", size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  labs(
    title = "Scatterplot of Residuals",
    x = "Residual EAC CCI",
    y = "Residual Current Strength",
    subtitle = paste("Pearson correlation:", round(residual_correlation2, 3), "| p-value:", signif(p_val3, 3))
  ) 


# Strength against climatology --------------------------------------------

time_range2 <- cop_model_data %>%
  reframe(time_range = range(date)) %>%
  pull()


clim_points2 <- tibble(date = seq(time_range2[1], time_range2[2], by = "1 day")) %>%
  mutate(doy = pmin(yday(date), 365)) %>%
  left_join(climatology_str2 %>% select(doy, str_clim), by = "doy")

clim_points2 <- clim_points2 %>% 
  mutate(date = as.POSIXct(date))





strength_with_clim2 <- cop_model_data %>%
  mutate(doy = yday(date)) %>%
  left_join(climatology_str2 %>% select(doy, str_clim), by = "doy") %>%
  mutate(
    vel_anomaly = mean_vel - str_clim,
    anom_label = case_when(
      vel_anomaly > 0 ~ "Anomaly (+)",
      TRUE ~ "Anomaly (-)")
  )


strength_with_clim2 <- strength_with_clim2 %>%
  mutate(date = as.POSIXct(date))


str_anom2 <- strength_with_clim2 %>% 
  select(vel_anomaly)

# Plot strength with climatology
p <- ggplot() +
  # Segments for anomalies
  geom_segment(
    data = strength_with_clim2,
    aes(x = date, xend = date, y = str_clim, yend = mean_vel,
        colour = anom_label, linetype = anom_label),
    linewidth = 1.2
  ) +
  # Climatology line
  geom_line(
    data = clim_points2,
    aes(x = date, y = str_clim, linetype = "Climatology", colour = "Climatology"),
    linewidth = 1.1
  ) +
  # Points for observed mean_vcur
  geom_point(
    data = strength_with_clim2,
    aes(x = date, y = mean_vel, fill = anom_label, colour = anom_label),
    size = 2.5, shape = 21
  ) +
  scale_colour_manual(
    breaks = c("Anomaly (+)", "Anomaly (-)", "Climatology"),
    values = c("Anomaly (+)" = "red", "Anomaly (-)" = "blue", "Climatology" = "black")
  ) +
  scale_fill_manual(
    breaks = c("Anomaly (+)", "Anomaly (-)"),
    values = c("Anomaly (+)" = "red", "Anomaly (-)" = "blue")
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
    title = "Current Strength Climatology (monthly average)"
  )

ggsave(file.path("output", "strength-with-clim2.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)
