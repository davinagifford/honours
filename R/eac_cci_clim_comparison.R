### eac_cci_clim_comparison.R
###
### plot the various climate indicators for comparison with the EAC CCI
###
### Created: 2025-08-14
### Author: Davina Gifford
### Last updated: 2025-08-14
### Edited by: Davina Gifford

library(tidyverse)
library(ggplot2)
library(lubridate)

# load SOI data

soi_data <- read_csv("data/soi_monthly.csv")

# Convert yyyymm to Date (first day of month)
soi_data <- soi_data %>%
  mutate(Date = ymd(paste0(Date, "01"))) %>% 
  filter(Date < "2025-01-01")

p <- ggplot(soi_data, aes(x = Date, y = SOI)) +
  geom_line(linewidth = 1, colour = "black") +
  geom_point(size = 3, shape = 21, fill = "black") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  geom_smooth(method = "loess", se = TRUE, color = "black", linetype = "dashed") +
  theme(axis.title = element_text(size = 14), 
        plot.title = element_text(size = 18, face = "bold")) +
  labs(x = "Time",
       y = "SOI",
       title = "Southern Oscillation Index ")

ggsave(file.path("output", "SOI.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)


# load SAM data

sam_data <- read_csv("data/AAO-SAM.csv")

sam_data <- sam_data %>%
  filter(year > 2009)

# Calculate monthly average
sam_monthly <- sam_data %>%
  group_by(year, month) %>%
  summarise(
    mean_sam = mean(aao_index_cdas, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    full_date = as.Date(paste(year, month, "01", sep = "-"))
  )

p <- ggplot(sam_monthly, aes(x = full_date, y = mean_sam)) +
  geom_line(linewidth = 1, colour = "black") +
  geom_point(size = 3, shape = 21, fill = "black") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  geom_smooth(method = "loess", se = TRUE, color = "black", linetype = "dashed") +
  theme(axis.title = element_text(size = 14), 
        plot.title = element_text(size = 18, face = "bold")) +
  labs(x = "Time",
       y = "SAM (Monthly Mean)",
       title = "Southern Annular Mode (SAM) - Monthly Average")

ggsave(file.path("output", "SAM.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)


# ENSO

 enso_data <- read_csv("data/nina34.anom.csv")

 enso_data <- enso_data %>% 
   filter(Date > "2010-01-01") %>% # Filter for dates after 2009
   filter(Date < "2025-01-01") %>% # Filter for dates before 2025
   mutate(anom = `Nino Anom 3.4 Index  using ersstv5 from CPC  missing value -99.99 https://psl.noaa.gov/data/timeseries/month/`) 
 
 head(enso_data)
 p <- ggplot(enso_data, aes(x = Date, y = anom )) +
   geom_line(linewidth = 1, colour = "black") +
   geom_point(size = 3, shape = 21, fill = "black") +
   scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
   geom_smooth(method = "loess", se = TRUE, color = "black", linetype = "dashed") +
   theme(axis.title = element_text(size = 14), 
         plot.title = element_text(size = 18, face = "bold")) +
   labs(x = "Time",
        y = "ENSO Index",
        title = "ENSO Index (NINO3.4 Anomaly)") 
 ggsave(file.path("output", "ENSO.png"),
        plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
        device = png)
 
 


# compare CCI values with SOI ---------------------------------------------

 

 
 # Join datasets by date
 combined_data <- month_data_t %>%
   mutate(Date = as.Date(trip_month)) %>%
   inner_join(soi_data, by = "Date")
 
 # Rescale SOI for plotting (adjust factor as needed)
 soi_scale_factor <- max(combined_data$eac_cci, na.rm = TRUE) / max(abs(combined_data$SOI), na.rm = TRUE)
 combined_data <- combined_data %>%
   mutate(SOI_scaled = SOI * soi_scale_factor)
 
 # Plot with dual y-axis
 p <- ggplot(combined_data, aes(x = Date)) +
   geom_line(aes(y = eac_cci), colour = "black", linewidth = 1) +
   geom_point(aes(y = eac_cci), size = 3, shape = 21, fill = "black") +
   geom_smooth(aes(y = eac_cci), method = "loess", se = TRUE, color = "black", linetype = "dashed") +
   geom_line(aes(y = SOI_scaled), colour = "blue", linewidth = 1) +
   geom_point(aes(y = SOI_scaled), size = 3, shape = 21, fill = "blue") +
   geom_smooth(aes(y = SOI_scaled), method = "loess", se = TRUE, color = "blue", linetype = "dashed") +
   scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
   scale_y_continuous(
     name = "EAC copepod composition index",
     sec.axis = sec_axis(~ . / soi_scale_factor, name = "SOI")
   ) +
   theme(
     axis.title.y = element_text(size = 14, color = "black"),
     axis.title.y.right = element_text(size = 14, color = "blue"),
     plot.title = element_text(size = 18, face = "bold")
   ) +
   labs(
     x = "Time",
     title = "EAC copepod composition index & SOI (monthly average)"
   )
 

 
 ggsave(file.path("output", "eac_cci_soi_combined.png"),
        plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
        device = png)

 
 # test similarity between cci and soi values
 
 
combined_data_forcorr <- combined_data %>%
   select(Date, eac_cci, SOI) %>%
   drop_na()
 
correlation <- cor(combined_data_forcorr$eac_cci, combined_data_forcorr$SOI, method = "pearson")
print(paste("Pearson correlation between EAC CCI and SOI:", correlation)) 

# anova
 
anova_result <- aov(eac_cci ~ SOI, data = combined_data_forcorr)

summary(anova_result)




# linear regression of EAC CCI against SOI

# Fit a linear model
model <- lm(eac_cci ~ SOI, data = combined_data)

# Summary of the model
summary(model)

# Plot the regression
plot(combined_data$SOI, combined_data$eac_cci, main = "Regression of EAC CCI on SOI",
     xlab = "SOI", ylab = "EAC CCI", pch = 19)
abline(model, col = "blue", lwd = 2)


plot(model$fitted.values, model$residuals,
     main = "Residuals vs Fitted Values",
     xlab = "Fitted Values",
     ylab = "Residuals",
     pch = 19, col = "darkgreen")
abline(h = 0, col = "red", lwd = 2)


# Save diagnostic plots to a PNG file
png(filename = "output/eac-soi-diagnostic.png", width = 1200, height = 600)

# Set up 2x2 layout and plot diagnostics
par(mfrow = c(2, 2))
plot(model)

# Close the PNG device
dev.off()


# Cross-correlation function to explore lagged relationships
ccf_result <- ccf(combined_data$SOI, combined_data$eac_cci, lag.max = 12, plot = TRUE,
                  main = "Cross-Correlation between SOI and EAC CCI")



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
combined_data$SOI_lagged <- dplyr::lag(combined_data$SOI, n = abs(best_lag))

# If lag is negative, shift EAC_CCI instead
if (best_lag < 0) {
  combined_data$eac_cci_lagged <- dplyr::lag(combined_data$eac_cci, n = abs(best_lag))
  lag_model <- lm(eac_cci_lagged ~ SOI, data = combined_data)
} else {
  lag_model <- lm(eac_cci ~ SOI_lagged, data = combined_data)
}

# View regression summary
summary(lag_model)



# Compare Index with Nina3.4 ----------------------------------------------

# Join datasets by date
combined_data_enso <- month_data_t %>%
  mutate(Date = as.Date(trip_month)) %>%
  inner_join(enso_data, by = "Date")

# Rescale for plotting (adjust factor as needed)
enso_scale_factor <- max(combined_data_enso$eac_cci, na.rm = TRUE) / max(abs(combined_data_enso$anom), na.rm = TRUE)
combined_data_enso <- combined_data_enso %>%
  mutate(enso_scaled = anom * enso_scale_factor)

# Plot with dual y-axis
p <- ggplot(combined_data_enso, aes(x = Date)) +
  geom_line(aes(y = eac_cci), colour = "black", linewidth = 1) +
  geom_point(aes(y = eac_cci), size = 3, shape = 21, fill = "black") +
  geom_smooth(aes(y = eac_cci), method = "loess", se = TRUE, color = "black", linetype = "dashed") +
  geom_line(aes(y = enso_scaled), colour = "blue", linewidth = 1) +
  geom_point(aes(y = enso_scaled), size = 3, shape = 21, fill = "blue") +
  geom_smooth(aes(y = enso_scaled), method = "loess", se = TRUE, color = "blue", linetype = "dashed") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(
    name = "EAC copepod composition index",
    sec.axis = sec_axis(~ . / enso_scale_factor, name = "ENSO Index (NINO3.4 Anomaly)")
  ) +
  theme(
    axis.title.y = element_text(size = 14, color = "black"),
    axis.title.y.right = element_text(size = 14, color = "blue"),
    plot.title = element_text(size = 18, face = "bold")
  ) +
  labs(
    x = "Time",
    title = "EAC copepod composition index & ENSO Index (NINO3.4 Anomaly)"
  )



ggsave(file.path("output", "eac_cci_enso_combined.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)


# test similarity between cci and enso values


enso_combined_data_forcorr <- combined_data_enso %>%
  select(Date, eac_cci, anom) %>%
  drop_na()

correlation_enso <- cor(enso_combined_data_forcorr$eac_cci, enso_combined_data_forcorr$anom, method = "pearson")
print(paste("Pearson correlation between EAC CCI and ENSO Index (Nino3.4 Anomaly:", correlation)) 

# anova

anova_result_enso <- aov(eac_cci ~ anom, data = enso_ombined_data_forcorr)

summary(anova_result_enso)




# linear regression of EAC CCI against SOI

# Fit a linear model
model_enso <- lm(eac_cci ~ anom, data = combined_data_enso)

# Summary of the model
summary(model_enso)

# Plot the regression
plot(combined_data_enso$anom, combined_data_enso$eac_cci, main = "Regression of EAC CCI on Nino3.4 Anomaly",
     xlab = "Nino3.4 Anomaly", ylab = "EAC CCI", pch = 19)
abline(model, col = "blue", lwd = 2)


plot(model_enso$fitted.values, model_enso$residuals,
     main = "Residuals vs Fitted Values",
     xlab = "Fitted Values",
     ylab = "Residuals",
     pch = 19, col = "darkgreen")
abline(h = 0, col = "red", lwd = 2)


# Save diagnostic plots to a PNG file
png(filename = "output/eac-enso-diagnostic.png", width = 1200, height = 600)

# Set up 2x2 layout and plot diagnostics
par(mfrow = c(2, 2))
plot(model_enso)

# Close the PNG device
dev.off()


# Cross-correlation function to explore lagged relationships
ccf_result_enso <- ccf(combined_data_enso$anom, combined_data_enso$eac_cci, lag.max = 12, plot = TRUE,
                  main = "Cross-Correlation between Nino3.4 Anomaly and EAC CCI")



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




# Create lagged ENSO variable
combined_data_enso$enso_lagged <- dplyr::lag(combined_data_enso$anom, n = abs(best_lag))

# If lag is negative, shift EAC_CCI instead
if (best_lag < 0) {
  combined_data_enso$eac_cci_lagged <- dplyr::lag(combined_data_enso$eac_cci, n = abs(best_lag))
  lag_model_enso <- lm(eac_cci_lagged ~ anom, data = combined_data_enso)
} else {
  lag_model_enso <- lm(eac_cci ~ enso_lagged, data = combined_data_enso)
}

# View regression summary
summary(lag_model_enso)



# Compare index with SAM --------------------------------------------------

# Join datasets by date
combined_data_sam <- month_data_t %>%
  mutate(full_date = as.Date(trip_month)) %>%
  inner_join(sam_data, by = "full_date")

# Rescale for plotting (adjust factor as needed)
sam_scale_factor <- max(combined_data_sam$eac_cci, na.rm = TRUE) / max(abs(combined_data_sam$aao_index_cdas), na.rm = TRUE)
combined_data_sam <- combined_data_sam %>%
  mutate(sam_scaled = aao_index_cdas * sam_scale_factor)

# Plot with dual y-axis
p <- ggplot(combined_data_sam, aes(x = full_date)) +
  geom_line(aes(y = eac_cci), colour = "black", linewidth = 1) +
  geom_point(aes(y = eac_cci), size = 3, shape = 21, fill = "black") +
  geom_smooth(aes(y = eac_cci), method = "loess", se = TRUE, color = "black", linetype = "dashed") +
  geom_line(aes(y = sam_scaled), colour = "blue", linewidth = 1) +
  geom_point(aes(y = sam_scaled), size = 3, shape = 21, fill = "blue") +
  geom_smooth(aes(y = sam_scaled), method = "loess", se = TRUE, color = "blue", linetype = "dashed") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(
    name = "EAC copepod composition index",
    sec.axis = sec_axis(~ . / sam_scale_factor, name = "Southern Annular Mode")
  ) +
  theme(
    axis.title.y = element_text(size = 14, color = "black"),
    axis.title.y.right = element_text(size = 14, color = "blue"),
    plot.title = element_text(size = 18, face = "bold")
  ) +
  labs(
    x = "Time",
    title = "EAC copepod composition index & Southern Annular Mode"
  )



ggsave(file.path("output", "eac_cci_sam_combined.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)


# test similarity between cci and enso values


sam_combined_data_forcorr <- combined_data_sam %>%
  select(full_date, eac_cci, aao_index_cdas) %>%
  drop_na()

correlation_sam <- cor(sam_combined_data_forcorr$eac_cci, sam_combined_data_forcorr$aao_index_cdas, method = "pearson")
print(paste("Pearson correlation between EAC CCI and Southern Annular Mode:", correlation)) 

# anova

anova_result_sam <- aov(eac_cci ~ aao_index_cdas, data = sam_combined_data_forcorr)

summary(anova_result_sam)




# linear regression of EAC CCI against SOI

# Fit a linear model
model_sam <- lm(eac_cci ~ aao_index_cdas, data = combined_data_sam)

# Summary of the model
summary(model_sam)

# Plot the regression
plot(combined_data_sam$aao_index_cdas, combined_data_sam$eac_cci, main = "Regression of EAC CCI on SAM",
     xlab = "SAM", ylab = "EAC CCI", pch = 19)
abline(model, col = "blue", lwd = 2)


plot(model_sam$fitted.values, model_sam$residuals,
     main = "Residuals vs Fitted Values",
     xlab = "Fitted Values",
     ylab = "Residuals",
     pch = 19, col = "darkgreen")
abline(h = 0, col = "red", lwd = 2)


# Save diagnostic plots to a PNG file
png(filename = "output/eac-sam-diagnostic.png", width = 1200, height = 600)

# Set up 2x2 layout and plot diagnostics
par(mfrow = c(2, 2))
plot(model_sam)

# Close the PNG device
dev.off()


# Cross-correlation function to explore lagged relationships
ccf_result <- ccf(combined_data_sam$aao_index_cdas, combined_data_sam$eac_cci, lag.max = 12, plot = TRUE,
                  main = "Cross-Correlation between SAM and EAC CCI")



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




# Create lagged SAM variable
combined_data_sam$sam_lagged <- dplyr::lag(combined_data_sam$aao_index_cdas, n = abs(best_lag))

# If lag is negative, shift EAC_CCI instead
if (best_lag < 0) {
  combined_data_sam$eac_cci_lagged <- dplyr::lag(combined_data_sam$eac_cci, n = abs(best_lag))
  lag_model_sam <- lm(eac_cci_lagged ~ aao_index_cdas, data = combined_data_sam)
} else {
  lag_model_sam <- lm(eac_cci ~ sam_lagged, data = combined_data_sam)
}

# View regression summary
summary(lag_model_sam)



