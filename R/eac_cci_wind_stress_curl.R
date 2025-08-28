### eac_cci_wind_stress_curl.R
###
### plot the wind stress curl for comparison with the EAC CCI
###
### Created: 2025-08-22
### Author: Davina Gifford
### Last updated: 2025-08-22
### Edited by: Davina Gifford

library(tidyverse)
library(ggplot2)
library(lubridate)
library(tidync)

# load wind stress data
# sourced from copernicus 

wind_stress <- "D:/HONOURS/DavinaG_2025_Honours/data/cmems_obs-wind_glo_phy_my_l4_P1M_1755838877871.nc"
if (!tolower(Sys.info()[["sysname"]]) == "sunos") {
  tnc <- tidync(wind_stress)
  print(tnc)
}

# explore variables
variables <- hyper_vars(tnc)
print(variables)


# extract eastward stress data
east_stress <- tnc %>%
  hyper_tibble(select_var = "eastward_stress")
print(east_stress)

# extract northward stress data
north_stress <- tnc %>%
  hyper_tibble(select_var = "northward_stress")
print(north_stress)

# Convert time to date 
east_stress <- east_stress %>%
  mutate(date = as.Date(time)) %>% 
  mutate(year = year(date), month = month(date)) %>% 
  group_by(year, month) %>%
  summarise(mean_stress = mean(eastward_stress, na.rm = TRUE), .groups = "drop")

north_stress <- north_stress %>% 
  mutate(date = as.Date(time)) %>% 
  mutate(year = year(date), month = month(date)) %>% 
  group_by(year, month) %>%
  summarise(mean_stress = mean(northward_stress, na.rm = TRUE), .groups = "drop")


# Combine eastward and northward stress into a single data frame

wind_stress_data <- east_stress %>%
  inner_join(north_stress, by = c("time", "longitude", "latitude")) %>% #join by time lat and long.
  mutate(date = as.Date(time)) %>% 
  mutate(month = floor_date(date, "month")) # add a month column



# Function to compute curl for one month
compute_curl <- function(data) {
  # Pivot to wide format for matrix operations
  tau_x <- data %>%
    select(latitude, longitude, eastward_stress) %>%
    pivot_wider(names_from = longitude, values_from = eastward_stress) %>%
    select(-latitude) %>%
    as.matrix()

  
  tau_y <- data %>%
    select(latitude, longitude, northward_stress) %>%
    pivot_wider(names_from = longitude, values_from = northward_stress) %>%
    select(-latitude) %>%
    as.matrix()
  
  # Compute gradients
  d_tau_y_dx <- apply(tau_y, 1, diff)  # across longitude
  d_tau_x_dy <- apply(tau_x, 2, diff)  # across latitude
  
  # Match dimensions by trimming
  min_rows <- min(nrow(d_tau_y_dx), nrow(d_tau_x_dy))
  min_cols <- min(ncol(d_tau_y_dx), ncol(d_tau_x_dy))
  
  curl <- d_tau_y_dx[1:min_rows, 1:min_cols] - d_tau_x_dy[1:min_rows, 1:min_cols]
  
  return(mean(curl, na.rm = TRUE))
}

# Apply by month
monthly_curl <- wind_stress_data %>%
  group_by(month) %>%
  group_modify(~ tibble(mean_curl = compute_curl(.x))) %>%
  ungroup()




p <- ggplot(monthly_curl, aes(x = month, y = mean_curl)) +
  geom_line(linewidth = 1, colour = "black") +
  geom_point(size = 3, shape = 21, fill = "black") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  geom_smooth(method = "loess", se = TRUE, color = "black", linetype = "dashed") +
  theme(axis.title = element_text(size = 14), 
        plot.title = element_text(size = 18, face = "bold")) +
  labs(x = "Time",
       y = "Mean Wind Stress Curl",
       title = "Wind stress curl - monthly average ")

ggsave(file.path("output", "mean-wind-stress.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)





# compare CCI values with SOI ---------------------------------------------

 

 
 # Join datasets by date
 combined_data_curl <- month_data_t %>%
   mutate(month = as.Date(trip_month)) %>%
   inner_join(monthly_curl, by = "month")
 
 # Rescale SOI for plotting (adjust factor as needed)
 curl_scale_factor <- max(combined_data_curl$eac_cci, na.rm = TRUE) / max(abs(combined_data_curl$mean_curl), na.rm = TRUE)
 combined_data_curl <- combined_data_curl %>%
   mutate(curl_scaled = mean_curl * curl_scale_factor)
 
 # Plot with dual y-axis
 p <- ggplot(combined_data_curl, aes(x = month)) +
   geom_line(aes(y = eac_cci), colour = "black", linewidth = 1) +
   geom_point(aes(y = eac_cci), size = 3, shape = 21, fill = "black") +
   geom_smooth(aes(y = eac_cci), method = "loess", se = TRUE, color = "black", linetype = "dashed") +
   geom_line(aes(y = curl_scaled), colour = "blue", linewidth = 1) +
   geom_point(aes(y = curl_scaled), size = 3, shape = 21, fill = "blue") +
   geom_smooth(aes(y = curl_scaled), method = "loess", se = TRUE, color = "blue", linetype = "dashed") +
   scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
   scale_y_continuous(
     name = "EAC copepod composition index",
     sec.axis = sec_axis(~ . / curl_scale_factor, name = "Mean Wind Stress Curl")
   ) +
   theme(
     axis.title.y = element_text(size = 14, color = "black"),
     axis.title.y.right = element_text(size = 14, color = "blue"),
     plot.title = element_text(size = 18, face = "bold")
   ) +
   labs(
     x = "Time",
     title = "EAC copepod composition index & Mean Wind Stress Curl (monthly average)"
   )
 

 
 ggsave(file.path("output", "eac_cci_curl_combined.png"),
        plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
        device = png)

 
 # test similarity between cci and wind curl values
 
 
curl_combined_data_forcorr <- combined_data_curl %>%
   select(month, eac_cci, mean_curl) %>%
   drop_na()
 
correlation_curl <- cor(curl_combined_data_forcorr$eac_cci, curl_combined_data_forcorr$mean_curl, method = "pearson")
print(paste("Pearson correlation between EAC CCI and SOI:", correlation_curl)) 

# anova
 
anova_result_curl <- aov(eac_cci ~ mean_curl, data = curl_combined_data_forcorr)

summary(anova_result_curl)




# linear regression of EAC CCI against SOI

# Fit a linear model
model_curl <- lm(eac_cci ~ mean_curl, data = combined_data_curl)

# Summary of the model
summary(model_curl)

# Plot the regression
plot(combined_data_curl$mean_curl, combined_data_curl$eac_cci, main = "Regression of EAC CCI on Mean wind stress curl",
     xlab = "Mean wind stress curl", ylab = "EAC CCI", pch = 19)
abline(model_curl, col = "blue", lwd = 2)


plot(model_curl$fitted.values, model_curl$residuals,
     main = "Residuals vs Fitted Values",
     xlab = "Fitted Values",
     ylab = "Residuals",
     pch = 19, col = "darkgreen")
abline(h = 0, col = "red", lwd = 2)


# Save diagnostic plots to a PNG file
png(filename = "output/eac-curl-diagnostic.png", width = 1200, height = 600)

# Set up 2x2 layout and plot diagnostics
par(mfrow = c(2, 2))
plot(model_curl)

# Close the PNG device
dev.off()


# Cross-correlation function to explore lagged relationships
ccf_result <- ccf(combined_data_curl$mean_curl, combined_data_curl$eac_cci, lag.max = 12, plot = TRUE,
                  main = "Cross-Correlation between Mean wind stress curl and EAC CCI")



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




# Create lagged SOI variable
combined_data_curl$curl_lagged <- dplyr::lag(combined_data_curl$mean_curl, n = abs(best_lag))

# If lag is negative, shift EAC_CCI instead
if (best_lag < 0) {
  combined_data_curl$eac_cci_lagged <- dplyr::lag(combined_data_curl$eac_cci, n = abs(best_lag))
  lag_model_curl <- lm(eac_cci_lagged ~ mean_curl, data = combined_data_curl)
} else {
  lag_model_curl <- lm(eac_cci ~ curl_lagged, data = combined_data_curl)
}

# View regression summary
summary(lag_model_curl)

