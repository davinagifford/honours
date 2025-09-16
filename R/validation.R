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

# read mooring data


# in tidy 

datafile <- "D:/HONOURS/DavinaG_2025_Honours/data/EAC_filled-daily-distance-depth-gridded-product_20120401-20220727.nc"
if (!tolower(Sys.info()[["sysname"]]) == "sunos") {
  tnc <- tidync(datafile)
  print(tnc)
}

# explore variables
variables <- hyper_vars(tnc)
print(variables)

# extract velocity data
data_tbl <- tnc %>%
  hyper_tibble(select_var = "VCUR")
print(data_tbl)


raw_nc_data <- tnc %>%
  hyper_tibble(select_var = "VCUR")

# look at the data before doing any mutates
raw_nc_data <- raw_nc_data %>% 
  filter(DEPTH < 1000) %>%
  filter(LONGITUDE > 153.9) %>% 
  filter(LONGITUDE < 155.0)

# filter data to desired depth and longitude range. Create monthly average velocity data
data_tbl <- data_tbl %>%
  filter(DEPTH < 1000) %>%
  filter(LONGITUDE > 153.9) %>% 
  filter(LONGITUDE < 155.0) %>% 
  mutate(
    date = as.Date(sub("T.*", "", TIME)),
    year = year(date),
    month = month(date)
  ) %>%
  group_by(year, month) %>%
  summarise(mean_vcur = mean(VCUR, na.rm = TRUE), .groups = "drop") %>% 
  mutate(date = as.Date(paste(year, month, "01", sep = "-"))) 

# convert velocity to positive strength value
data_tbl <- data_tbl %>%
  mutate(mean_vcur = abs(mean_vcur))


# plot the data
p <- ggplot(data_tbl, aes(x = date, y = mean_vcur)) +
  geom_line(linewidth = 1, colour = "black") +
  geom_point(size = 3, shape = 21, fill = "black") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  geom_smooth(method = "loess", se = TRUE, color = "black", linetype = "dashed") +
  labs(
    x = "Time",
    y = expression("Monthly Mean strength (m/s"^{-1}*")"),
    title = "Mean Strength Over Time (100m depth)"
  ) +
  theme(axis.title = element_text(size = 14), 
        plot.title = element_text(size = 18, face = "bold"))

ggsave(file.path("output", "mean-velocity_100m.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)


# extract temperature data
temp_tbl <- tnc %>%
  hyper_tibble(select_var = "TEMP")
print(temp_tbl)

# filter data to desired depth and longitude range. Create monthly average velocity data
temp_tbl <- temp_tbl %>%
  filter(DEPTH < 1000) %>%
  filter(LONGITUDE > 153.9) %>% 
  filter(LONGITUDE < 155.0) %>% 
  mutate(
    date = as.Date(sub("T.*", "", TIME)),
    year = year(date),
    month = month(date)
  ) %>%
  group_by(year, month) %>%
  summarise(mean_temp = mean(TEMP, na.rm = TRUE), .groups = "drop") %>% 
  mutate(date = as.Date(paste(year, month, "01", sep = "-"))) 



# plot the data
p <- ggplot(temp_tbl, aes(x = date, y = mean_temp)) +
  geom_line(linewidth = 1, colour = "black") +
  geom_point(size = 3, shape = 21, fill = "black") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  geom_smooth(method = "loess", se = TRUE, color = "black", linetype = "dashed") +
  labs(
    x = "Time",
    y = expression("Monthly Mean Temperature"),
    title = "Mean Temperature Over Time (100m depth)"
  ) +
  theme(axis.title = element_text(size = 14), 
        plot.title = element_text(size = 18, face = "bold"))

ggsave(file.path("output", "mean-temp_100m.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)

# calculate a climatology


current_data <- data_tbl %>%
  mutate(
    doy = yday(date)
  )

lm_fit_strength <- gam(mean_vcur ~ s(doy, bs = "cc", k = 5),
                       data = current_data,
                       method = "REML",
                       knots = list(doy = c(0, 365)))

mean_time <- median(current_data$date)
time0 <- min(current_data$date)

climatology_str <-
  tibble(date = floor_date(mean_time, unit = "year") + days(0:364),
         doy = yday(date),
         time_x = time_length(interval(time0, mean_time), unit = "day"))

clim_terms_str <- predict(lm_fit_strength, type = "terms", newdata = climatology_str)

head(clim_terms_str)

climatology_str <-
  climatology_str %>%
  mutate(doy_eff = clim_terms_str[, "s(doy)"],
         intercept = coef(lm_fit_strength)["(Intercept)"],
         str_clim = doy_eff + intercept) %>%
  select(!c(doy_eff, intercept))

climatology_str <- climatology_str %>%
  mutate(date = as.POSIXct(date))

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

# compare monthly climatology values---------------


# plot strength climatology with cci climatology
climatology_str <- climatology_str %>%
  mutate(clim_date = as.Date(date))

climatology <- climatology %>%
  mutate(clim_date = as.Date(sample_time))

month_climatology <- climatology %>% 
  mutate(month = month(clim_date)) %>%
  group_by(month) %>%
  summarise(month_clim = mean(eac_cci, na.rm = TRUE), .groups = "drop")

month_clim_str <- climatology_str %>%
  mutate(month = month(clim_date)) %>%
  group_by(month) %>%
  summarise(month_clim_str = mean(str_clim, na.rm = TRUE), .groups = "drop")

clim_compare_short <- month_clim_str %>%
  select(month, month_clim_str) %>%
  left_join(
    month_climatology,
    by = "month"
  )



# Pivot longer for easier plotting
clim_compare_long <- clim_compare_short %>%
  pivot_longer(cols = c(month_clim_str, month_clim),
               names_to = "climatology_type",
               values_to = "value")

# Calculate correlation


# Run Pearson correlation
correlation <- cor(clim_compare_short$month_clim, clim_compare_short$month_clim_str, method = "pearson", use = "complete.obs")
cor_test <- cor.test(clim_compare_short$month_clim, clim_compare_short$month_clim_str, method = "pearson", use = "complete.obs")

# Output results
print(paste("Pearson correlation:", round(correlation, 3)))
print(cor_test)

# run 


p_val <- cor_test$p.value

# run Kendall correlation
correlation_k <- cor(clim_compare_short$month_clim, clim_compare_short$month_clim_str, method = "kendall", use = "complete.obs")
cor_test_k <- cor.test(clim_compare_short$month_clim, clim_compare_short$month_clim_str, method = "kendall", use = "complete.obs")

# Output results
print(paste("Kendall correlation:", round(correlation_k, 3)))
print(cor_test_k)


p_val_k <- cor_test_k$p.value



# Plot comparison
p <- ggplot(clim_compare_long, aes(x = month, y = value, color = climatology_type)) +
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
    subtitle = paste("Kendall correlation:", round(correlation_k, 3), "| p-value:", signif(p_val_k, 3))
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )


ggsave(file.path("output", "climatology-comparison_kendall.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)


# plot a scatterplot of the two

ggplot(data = clim_compare_short) +
  geom_point(mapping = aes(x = month_clim, y = month_clim_str)) +
  geom_smooth(mapping = aes(x = month_clim, y = month_clim_str), method = "lm", se = TRUE, color = "blue", linetype = "dashed") +
  labs(
    x = "EAC CCI Climatology",
    y = "Current Strength Climatology",
    title = "Scatterplot of Current Strength vs EAC CCI Climatologies",
    subtitle = paste("Pearson correlation:", round(correlation, 3), "| p-value:", signif(p_val, 3))
  ) 


# correlate the residuals

# Fit a linear model
str_model <- lm(month_clim ~ month_clim_str, data = clim_compare_short)

# Summary of the model
summary(str_model)



residuals_str <- resid(str_model)
predicted_str <- fitted(str_model)



plot(clim_compare_short$month_clim_str, str_model$residuals,
     main = "Residuals vs Fitted Values",
     xlab = "Fitted Values",
     ylab = "Residuals",
     pch = 19, col = "darkgreen")
abline(h = 0, col = "red", lwd = 2)




# # is there a relationship between velocity and temperature?  ------------



vel_temp <- data_tbl %>%
  left_join(temp_tbl, by = "date") %>%
  select(date, mean_vcur, mean_temp)

# Rescale for plotting (adjust factor as needed)
temp_scale_factor <- max(vel_temp$mean_vcur, na.rm = TRUE) / max(abs(vel_temp$mean_temp), na.rm = TRUE)
vel_temp <- vel_temp %>%
  mutate(temp_scaled = mean_temp * temp_scale_factor)

# Plot with dual y-axis
p <- ggplot(vel_temp, aes(x = date)) +
  geom_line(aes(y = mean_vcur), colour = "black", linewidth = 1) +
  geom_point(aes(y = mean_vcur), size = 3, shape = 21, fill = "black") +
  geom_smooth(aes(y = mean_vcur), method = "loess", se = TRUE, color = "black", linetype = "dashed") +
  geom_line(aes(y = temp_scaled), colour = "blue", linewidth = 1) +
  geom_point(aes(y = temp_scaled), size = 3, shape = 21, fill = "blue") +
  geom_smooth(aes(y = temp_scaled), method = "loess", se = TRUE, color = "blue", linetype = "dashed") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(
    name = "Mean Monthly strength",
    sec.axis = sec_axis(~ . / temp_scale_factor, name = "Mean monthly Temperature")
  ) +
  theme(
    axis.title.y = element_text(size = 14, color = "black"),
    axis.title.y.right = element_text(size = 14, color = "blue"),
    plot.title = element_text(size = 18, face = "bold")
  ) +
  labs(
    x = "Time",
    title = "Mean monthly temperature and mean monthly strength (100m depth)"
  )



ggsave(file.path("output", "vel_temp_combined.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)


# test similarity between cci and enso values


vt_combined_data_forcorr <- vel_temp %>%
  select(date, mean_vcur, mean_temp) %>%
  drop_na()

correlation_vt <- cor(vt_combined_data_forcorr$mean_vcur, vt_combined_data_forcorr$mean_temp, method = "pearson")
print(paste("Pearson correlation between Mean Temp and Mean Velocity:", correlation_vt))

# anova

anova_result_vt <- aov(mean_vcur ~ mean_temp, data = vel_temp)

summary(anova_result_vt)




# linear regression of mean strength against mean temperature

# Fit a linear model
model_vt <- lm(mean_vcur ~ mean_temp, data = vel_temp)

# Summary of the model
summary(model_vt)



# Save diagnostic plots to a PNG file
png(filename = "output/vel-temp-diagnostic.png", width = 1200, height = 600)

# Set up 2x2 layout and plot diagnostics
par(mfrow = c(2, 2))
plot(model_vt)

# Close the PNG device
dev.off()


# Cross-correlation function to explore lagged relationships
ccf_result <- ccf(vel_temp$mean_vcur, vel_temp$mean_temp, lag.max = 12, plot = TRUE,
                  main = "Cross-Correlation between Mean Monthly velocity and Mean Monthly Temperature")



# Extract correlation values and corresponding lags
correlations <- ccf_result$acf[,1,1]
lags <- ccf_result$lag[,1,1]

# Find the lag with the highest absolute correlation
max_index <- which.max(abs(correlations))
best_lag <- lags[max_index]
best_corr <- correlations[max_index]

# Print result
cat("Lag with highest cross-correlation:", best_lag, "months\n")
cat("Correlation coefficient:", round(best_corr, 3), "\n")




# Create lagged temp variable
vel_temp$temp_lagged <- dplyr::lag(vel_temp$mean_temp, n = abs(best_lag))

# If lag is negative, shift strength instead
if (best_lag < 0) {
  vel_temp$vel_lagged <- dplyr::lag(vel_temp$mean_vcur, n = abs(best_lag))
  lag_model_vt <- lm(vel_lagged ~ mean_temp, data = vel_temp)
} else {
  lag_model_vt <- lm(mean_vcur ~ temp_lagged, data = vel_temp)
}

# View regression summary
summary(lag_model_vt)





# Compare strength with index ---------------------------------------------



# plot with the index values

with_vel <- anomaly %>%
  mutate(date = as.Date(trip_month)) %>%
  full_join(data_tbl, by = "date")

velocity_scale <- max(with_vel$anomaly, na.rm = TRUE) / max(abs(with_vel$mean_vcur), na.rm = TRUE)

with_vel <- with_vel %>%
  mutate(vel_scale = mean_vcur * velocity_scale)

p <- ggplot(with_vel, aes(x = date)) +
  geom_line(aes(y = anomaly), colour = "black", linewidth = 1) +
  geom_point(aes(y = anomaly), size = 3, shape = 21, fill = "black") +
  geom_smooth(aes(y = anomaly), method = "loess", se = TRUE, color = "black", linetype = "dashed") +
  geom_line(aes(y = vel_scale), colour = "blue", linewidth = 1) +
  geom_point(aes(y = vel_scale), size = 3, shape = 21, fill = "blue") +
  geom_smooth(aes(y = vel_scale), method = "loess", se = TRUE, color = "blue", linetype = "dashed") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(
    name = "EAC copepod composition index - anomaly",
    sec.axis = sec_axis(~ . / velocity_scale, name = "Monthly mean strength")
  ) +
  theme(
    axis.title.y = element_text(size = 14, color = "black"),
    axis.title.y.right = element_text(size = 14, color = "blue"),
    plot.title = element_text(size = 18, face = "bold")
  ) +
  labs(
    x = "Time",
    title = "EAC copepod composition index & Monthly mean strength"
  )

ggsave(file.path("output", "eac_cci_vel_combined.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = "png")


# test similarity between cci and strength values

with_vel <- with_vel %>% 
  mutate(
    year = year(trip_month),
    month = month(trip_month)
  ) 

str_combined_data_forcorr <- with_vel %>%
  select(month, observed_eac_cci, mean_vcur) %>%
  drop_na()

correlation_str <- cor(str_combined_data_forcorr$observed_eac_cci, str_combined_data_forcorr$mean_vcur, method = "pearson")
print(paste("Pearson correlation between EAC CCI and strength:", correlation_str)) 

# anova

anova_result_str <- aov(observed_eac_cci ~ mean_vcur, data = str_combined_data_forcorr)

summary(anova_result_str)

# linear regression of EAC CCI against strength

# Fit a linear model
model_str <- lm(observed_eac_cci ~ mean_vcur, data = with_vel)

# Summary of the model
summary(model_str)

model_res = resid(model_str)

anova(model_str)



ggplot(with_vel, aes(x = mean_vcur, y = observed_eac_cci)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(title = "Linear Regression of EAC CCI on Strength",
       x = expression("Monthly Mean strength (m/s"^{-1}*")"),
       y = "EAC CCI ") +
  theme_minimal()


# Plot the regression
plot(with_vel$mean_vcur, with_vel$observed_eac_cci, main = "Regression of EAC CCI on strength",
     xlab = "Mean Strength", ylab = "EAC CCI", pch = 19)
abline(model_str, col = "blue", lwd = 2)


plot(with_vel$observed_eac_cci, model_res,
     main = "Residuals vs EAC CCI",
     xlab = "EAC CCI",
     ylab = "Residuals",
     pch = 19, col = "darkgreen")
abline(h = 0, col = "red", lwd = 2)


# Save diagnostic plots to a PNG file
png(filename = "output/eac-str-diagnostic.png", width = 1200, height = 600)

# Set up 2x2 layout and plot diagnostics
par(mfrow = c(2, 2))
plot(model_str)

# Close the PNG device
dev.off()


# Remove rows with NA in either column
clean_data <- na.omit(with_vel[, c("mean_vcur", "observed_eac_cci")])

# Run CCF on cleaned data
ccf_result_str <- ccf(clean_data$mean_vcur, clean_data$observed_eac_cci, lag.max = 12, plot = TRUE,
                      main = "Cross-Correlation between Strength and EAC CCI")



# Extract correlation values and corresponding lags
correlations <- ccf_result_str$acf[,1,1]
lags <- ccf_result_str$lag[,1,1]

# Find the lag with the highest absolute correlation
max_index <- which.max(abs(correlations))
best_lag <- lags[max_index]
best_corr <- correlations[max_index]

# Print result
cat("Lag with highest cross-correlation:", best_lag, "months\n")
cat("Correlation coefficient:", round(best_corr, 3), "\n")




# Create lagged strength variable
with_vel$str_lagged <- dplyr::lag(with_vel$mean_vcur, n = abs(best_lag))

# If lag is negative, shift EAC_CCI instead
if (best_lag < 0) {
  with_vel$eac_cci_lagged <- dplyr::lag(with_vel$observed_eac_cci, n = abs(best_lag))
  lag_model_str <- lm(eac_cci_lagged ~ mean_vcur, data = with_vel)
} else {
  lag_model_str <- lm(observed_eac_cci ~ str_lagged, data = with_vel)
}

# View regression summary
summary(lag_model_str)


cor.test(with_vel$observed_eac_cci, with_vel$str_lagged, method = "pearson", use = "complete.obs")



# compare anomaly of index with anomaly of strength -----------------------

cci_strn <- anomaly %>%
  select(date = trip_month, cci_anom = anomaly) %>%
  left_join(str_anom %>% mutate(date = data_tbl$date), by = "date")

anom_combined_data_forcorr <- cci_strn %>%
  select(date, cci_anom, vcur_anomaly) %>%
  drop_na()

correlation_anom <- cor(anom_combined_data_forcorr$cci_anom, anom_combined_data_forcorr$vcur_anomaly, method = "pearson")
print(paste("Pearson correlation between EAC CCI and strength anomalies:", correlation_anom)) 

# anova

anova_result_anom <- aov(cci_anom ~ vcur_anomaly, data = anom_combined_data_forcorr)

summary(anova_result_anom)

# linear regression of EAC CCI against strength anomalies

# Fit a linear model
model_anom <- lm(cci_anom ~ vcur_anomaly, data = anom_combined_data_forcorr)

# Summary of the model
summary(model_anom)


ggplot(anom_combined_data_forcorr, aes(x = vcur_anomaly, y = cci_anom)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(title = "Linear Regression of EAC CCI anomalies on Strength aomalies",
       x = expression("Anomaly of Monthly Mean strength (m/s"^{-1}*")"),
       y = "EAC CCI Anomalies ") +
  theme_minimal()

# Plot the regression
plot(anom_combined_data_forcorr$vcur_anomaly, anom_combined_data_forcorr$cci_anom, main = "Regression of EAC CCI anomalies on strength anomalies",
     xlab = "Mean Strength anomalies", ylab = "EAC CCI anomalies", pch = 19)
abline(model_anom, col = "blue", lwd = 2)


plot(model_anom$fitted.values, model_anom$residuals,
     main = "Residuals vs Fitted Values",
     xlab = "Fitted Values",
     ylab = "Residuals",
     pch = 19, col = "darkgreen")
abline(h = 0, col = "red", lwd = 2)


# Save diagnostic plots to a PNG file
png(filename = "output/anom-diagnostic.png", width = 1200, height = 600)

# Set up 2x2 layout and plot diagnostics
par(mfrow = c(2, 2))
plot(model_anom)

# Close the PNG device
dev.off()


# Remove rows with NA in either column
clean_data <- na.omit(with_vel[, c("mean_vcur", "observed_eac_cci")])

# Run CCF on cleaned data
ccf_result_anom <- ccf(anom_combined_data_forcorr$vcur_anomaly, anom_combined_data_forcorr$cci_anom, lag.max = 12, plot = TRUE,
                       main = "Cross-Correlation between the anomalies in Strength and EAC CCI")



# Extract correlation values and corresponding lags
correlations <- ccf_result_anom$acf[,1,1]
lags <- ccf_result_anom$lag[,1,1]

# Find the lag with the highest absolute correlation
max_index <- which.max(abs(correlations))
best_lag <- lags[max_index]
best_corr <- correlations[max_index]

# Print result
cat("Lag with highest cross-correlation:", best_lag, "months\n")
cat("Correlation coefficient:", round(best_corr, 3), "\n")




# Create lagged strength variable
anom_combined_data_forcorr$vcur_lagged <- dplyr::lag(anom_combined_data_forcorr$vcur_anomaly, n = abs(best_lag))

# If lag is negative, shift EAC_CCI instead
if (best_lag < 0) {
  anom_combined_data_forcorr$cci_anom_lagged <- dplyr::lag(anom_combined_data_forcorr$cci_anom, n = abs(best_lag))
  lag_model_anom <- lm(cci_anom_lagged ~ vcur_anomaly, data = anom_combined_data_forcorr)
} else {
  lag_model_anom <- lm(cci_anom ~ vcur_lagged, data = anom_combined_data_forcorr)
}

# View regression summary
summary(lag_model_anom)


ggplot(anom_combined_data_forcorr, aes(x = vcur_anomaly, y = cci_anom)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(
    x = "EAC Strength Anomaly (VCUR)",
    y = "EAC CCI Anomaly",
    title = "Relationship Between EAC Strength and CCI Anomalies"
  )

# compare temp with sst ---------------------------------------------------

# read sst data
sst <- samples %>% 
  select(sample_time, sst)
colnames(sst) <- c("date", "sst")
sst <- sst %>%
  mutate(
    date = as.Date(sub("T.*", "", date)),
    year = year(date),
    month = month(date)
  ) %>%
  group_by(year, month) %>%
  summarise(sst = mean(sst, na.rm = TRUE), .groups = "drop") %>% 
  mutate(date = as.Date(paste(year, month, "01", sep = "-"))) 


sst_temp <- sst %>% 
  left_join(vel_temp %>% mutate(date = vel_temp$date), by = "date")

sst_temp <- sst_temp %>% 
  select(year, month, date, sst, mean_temp)


correlation_sst <- cor(sst_temp$sst, sst_temp$mean_temp, method = "pearson", use = "complete.obs")
print(paste("Pearson correlation between Index SST and observed Temperature:", correlation_sst)) 

# anova

anova_result_sst <- aov(sst ~ mean_temp, data = sst_temp)

summary(anova_result_sst)

# linear regression of EAC CCI against strength anomalies

# Fit a linear model
model_sst <- lm(sst ~ mean_temp, data = sst_temp)

# Summary of the model
summary(model_sst)

ggplot(sst_temp, aes(x = mean_temp, y = sst)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(title = "Linear Regression of SST on Mean monthly Temperatures",
       x = "Mean monthly Temperatures",
       y = "SST") +
  theme_minimal()


# Plot the regression
plot(anom_combined_data_forcorr$vcur_anomaly, anom_combined_data_forcorr$cci_anom, main = "Regression of EAC CCI anomalies on strength anomalies",
     xlab = "Mean Strength anomalies", ylab = "EAC CCI anomalies", pch = 19)
abline(model_anom, col = "blue", lwd = 2)


plot(model_anom$fitted.values, model_anom$residuals,
     main = "Residuals vs Fitted Values",
     xlab = "Fitted Values",
     ylab = "Residuals",
     pch = 19, col = "darkgreen")
abline(h = 0, col = "red", lwd = 2)


# Save diagnostic plots to a PNG file
png(filename = "output/anom-diagnostic.png", width = 1200, height = 600)

# Set up 2x2 layout and plot diagnostics
par(mfrow = c(2, 2))
plot(model_anom)

# Close the PNG device
dev.off()


# Remove rows with NA in either column
clean_data <- na.omit(with_vel[, c("mean_vcur", "observed_eac_cci")])

# Run CCF on cleaned data
ccf_result_anom <- ccf(anom_combined_data_forcorr$vcur_anomaly, anom_combined_data_forcorr$cci_anom, lag.max = 12, plot = TRUE,
                       main = "Cross-Correlation between the anomalies in Strength and EAC CCI")



# Extract correlation values and corresponding lags
correlations <- ccf_result_anom$acf[,1,1]
lags <- ccf_result_anom$lag[,1,1]

# Find the lag with the highest absolute correlation
max_index <- which.max(abs(correlations))
best_lag <- lags[max_index]
best_corr <- correlations[max_index]

# Print result
cat("Lag with highest cross-correlation:", best_lag, "months\n")
cat("Correlation coefficient:", round(best_corr, 3), "\n")




# Create lagged strength variable
anom_combined_data_forcorr$vcur_lagged <- dplyr::lag(anom_combined_data_forcorr$vcur_anomaly, n = abs(best_lag))

# If lag is negative, shift EAC_CCI instead
if (best_lag < 0) {
  anom_combined_data_forcorr$cci_anom_lagged <- dplyr::lag(anom_combined_data_forcorr$cci_anom, n = abs(best_lag))
  lag_model_anom <- lm(cci_anom_lagged ~ vcur_anomaly, data = anom_combined_data_forcorr)
} else {
  lag_model_anom <- lm(cci_anom ~ vcur_lagged, data = anom_combined_data_forcorr)
}

# View regression summary
summary(lag_model_anom)




# compare SST with Strength -----------------------------------------------

sst_str <- sst %>% 
  left_join(vel_temp %>% mutate(date = vel_temp$date), by = "date")

sst_str <- sst_str %>% 
  select(year, month, date, sst, mean_vcur)

correlation_sst_str <- cor(sst_str$sst, sst_str$mean_vcur, method = "pearson", use = "complete.obs")
print(paste("Pearson correlation between Index SST and observed strength:", correlation_sst_str)) 

# anova

anova_result_sst_str <- aov(sst ~ mean_vcur, data = sst_str)

summary(anova_result_sst_str)

# linear regression of EAC CCI against strength anomalies

# Fit a linear model
model_sst_str <- lm(sst ~ mean_vcur, data = sst_str)

# Summary of the model
summary(model_sst_str)


ggplot(sst_str, aes(x = mean_vcur, y = sst)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(title = "Linear Regression of SST on Mean monthly strength",
       x = expression("Anomaly of Monthly Mean strength (m/s"^{-1}*")"),
       y = "SST ") +
  theme_minimal()

# Remove rows with NA in either column
clean_data_sst_str <- na.omit(sst_str[, c("mean_vcur", "sst")])

# Run CCF on cleaned data
ccf_result_sst_str <- ccf(clean_data_sst_str$mean_vcur, clean_data_sst_str$sst, lag.max = 12, plot = TRUE,
                          main = "Cross-Correlation between the Mean strength and SST used in Index")



# Extract correlation values and corresponding lags
correlations <- ccf_result_sst_str$acf[,1,1]
lags <- ccf_result_sst_str$lag[,1,1]

# Find the lag with the highest absolute correlation
max_index <- which.max(abs(correlations))
best_lag <- lags[max_index]
best_corr <- correlations[max_index]

# Print result
cat("Lag with highest cross-correlation:", best_lag, "months\n")
cat("Correlation coefficient:", round(best_corr, 3), "\n")




# Create lagged strength variable
sst_str$vcur_lagged <- dplyr::lag(sst_str$mean_vcur, n = abs(best_lag))

# If lag is negative, shift SST instead
if (best_lag < 0) {
  sst_str$sst_lagged <- dplyr::lag(sst_str$sst, n = abs(best_lag))
  lag_model_sst_str <- lm(sst_lagged ~ mean_vcur, data = sst_str)
} else {
  lag_model_sst_str <- lm(sst ~ vcur_lagged, data = sst_str)
}

# View regression summary
summary(lag_model_sst_str)

cor(sst_str$sst_lagged, sst_str$mean_vcur, method = "pearson", use = "complete.obs")
cor.test(sst_str$sst_lagged, sst_str$mean_vcur, method = "pearson", use = "complete.obs")


# compare all the things --------------------------------------------------

# Combine temp, sst, strength, and index data

all_data <- with_vel %>%
  select(date, observed_eac_cci, mean_vcur) %>%
  left_join(vel_temp %>% select(date, mean_temp), by = "date") %>%
  left_join(sst %>% select(date, sst), by = "date")

# create a scatterplot matrix of the variables

pairs(all_data[, -1], main = "Scatterplot Matrix of EAC CCI, Strength, Temperature, and SST")
# calculate correlation matrix
cor_matrix <- cor(all_data[, -1], use = "pairwise.complete.obs")
print(cor_matrix)
# visualize correlation matrix
library(corrplot)
corrplot(cor_matrix, method = "circle", type = "upper", tl.col = "black", tl.srt = 45,
         title = "Correlation Matrix of EAC CCI, Strength, Temperature, and SST", mar = c(0,0,1,0))


cor(all_data$sst, all_data$observed_eac_cci, method = "pearson", use = "complete.obs")

summary(lm(sst ~ observed_eac_cci, data = all_data))
summary(lm(sst ~ mean_temp, data = all_data))
summary(lm(sst ~ mean_vcur, data = all_data))
summary(lm(mean_temp ~ observed_eac_cci, data = all_data))
summary(lm(mean_temp ~ mean_vcur, data = all_data))
summary(lm(mean_vcur ~ observed_eac_cci, data = all_data))
summary(lm(observed_eac_cci ~ mean_vcur, data = all_data))
summary(lm(observed_eac_cci ~ mean_temp, data = all_data))
summary(lm(observed_eac_cci ~ sst, data = all_data))
summary(lm(mean_vcur ~ sst, data = all_data))
summary(lm(mean_temp ~ sst, data = all_data))
summary(lm(mean_vcur ~ mean_temp, data = all_data))
summary(lm(mean_vcur ~ mean_temp, data = all_data))
summary(lm(north_index ~ mean_vcur, data = all_data))






# Validate the strength against just the index values from the North

all_data <- all_data %>% 
  left_join(month_data_north %>% select(date = trip_month, north_index = eac_cci), by = "date")


pairs(all_data[, -1], main = "Scatterplot Matrix of EAC CCI, Strength, Temperature, and SST")
# calculate correlation matrix
cor_matrix <- cor(all_data[, -1], use = "pairwise.complete.obs")
print(cor_matrix)


# Run Pearson correlation
correlation <- cor(all_data$mean_vcur, all_data$north_index, method = "pearson", use = "pairwise.complete.obs")
cor_test <- cor.test(all_data$mean_vcur, all_data$north_index, method = "pearson", use = "pairwise.complete.obs")

# Output results
print(paste("Pearson correlation:", round(correlation, 3)))
print(cor_test)

# Run kendall correlation
correlation_k <- cor(all_data$mean_vcur, all_data$north_index, method = "kendall", use = "pairwise.complete.obs")
cor_test_k <- cor.test(all_data$mean_vcur, all_data$north_index, method = "kendall", use = "pairwise.complete.obs")

# Output results
print(paste("Kendall correlation:", round(correlation_k, 3)))
print(cor_test_k)



# simple plot of Index against strength

# Remove rows with NA in either column
scatter_data <- with_vel %>%
  select(observed_eac_cci, mean_vcur) %>%
  drop_na()

# Calculate correlation
cor_result <- cor.test(scatter_data$observed_eac_cci, scatter_data$mean_vcur, method = "pearson")
corr_label <- sprintf("r = %.2f\np = %.3g", cor_result$estimate, cor_result$p.value)

# Plot
p <- ggplot(scatter_data, aes(x = mean_vcur, y = observed_eac_cci)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    x = "Observed Monthly Mean Current Strength (m/s)",
    y = "EAC Index (CCI)",
    title = "EAC Index vs Observed Current Strength"
  ) +
  annotate("text",
           x = min(scatter_data$mean_vcur, na.rm = TRUE),
           y = max(scatter_data$observed_eac_cci, na.rm = TRUE),
           label = corr_label,
           hjust = 0, vjust = 1, size = 5) +
  theme_minimal()

ggsave(file.path("output", "eac_vs_strength_scatter.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)


