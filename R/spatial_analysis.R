### spatial_analysis.R
###
### Analysis of the differences between north and south of the EAC.
###
### Created: 2025-07-04
### Author: Davina Gifford
### Last updated: 2025-08-30
### Edited by: Davina Gifford

# load libraries
library(dplyr)
library(lubridate)
library(ggplot2)

#Calculate monthly averages for south
monthly_avg_south <- month_data_south %>%
  mutate(month = month(trip_month, label = TRUE)) %>%  # Extract month as a factor with labels
  group_by(month) %>%
  summarise(avg_eac_cci = mean(eac_cci, na.rm = TRUE), .groups = "drop")

print(monthly_avg_south)

#Calculate monthly averages for north
monthly_avg_north <- month_data_north %>%
  mutate(month = month(trip_month, label = TRUE)) %>%  # Extract month as a factor with labels
  group_by(month) %>%
  summarise(avg_eac_cci = mean(eac_cci, na.rm = TRUE), .groups = "drop")

print(monthly_avg_north)


# Ensure both data frames are sorted by month
monthly_avg_south <- monthly_avg_south %>% arrange(month)
monthly_avg_north <- monthly_avg_north %>% arrange(month)

# Extract numeric vectors
south_values <- monthly_avg_south$avg_eac_cci
north_values <- monthly_avg_north$avg_eac_cci

# Run Pearson correlation
correlation <- cor(south_values, north_values, method = "pearson")
cor_test <- cor.test(south_values, north_values, method = "pearson")

# Output results
print(paste("Pearson correlation:", round(correlation, 3)))
print(cor_test)


p_val <- cor_test$p.value


# Plot the correlation

# make the dataframe

df <- data.frame(
  Month = monthly_avg_north$month,
  South = as.numeric(south_values),
  North = as.numeric(north_values)
)


comparison_ns_model <- lm(North ~ South, data = df)

# Save diagnostic plots to a PNG file
png(filename = "output/north-south-resid.png", width = 1200, height = 600)

# Set up 2x2 layout and plot diagnostics
par(mfrow = c(2, 2))
plot(comparison_ns_model)

# Close the PNG device
dev.off()




# Plot
ggplot(df, aes(x = South, y = North)) +
  geom_point(color = "black", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  theme(axis.title = element_text(size = 14), 
        plot.title = element_text(size = 18, face = "bold")) +
  labs(
    title = "Correlation of Monthly Average EAC CCI (South vs North)",
    x = "Southern Region EAC CCI",
    y = "Northern Region EAC CCI",
    subtitle = paste("Pearson correlation:", round(correlation, 3),
                     ", p-value = ", signif(p_val, 3))
  ) 
ggsave("output/eac_cci_correlation_north_south.png", width = 8, height = 6, dpi = 300)


# Assuming south_data and north_data have matching dates or trips
raw_df <- inner_join(south_data, north_data, by = "trip_month", suffix = c("_south", "_north"))
correlation_raw <- cor(raw_df$eac_cci_south, raw_df$eac_cci_north, method = "pearson")
cor_test_raw <- cor.test(raw_df$eac_cci_south, raw_df$eac_cci_north, method = "pearson")
print(correlation_raw)
print(cor_test_raw)


# Autocorrelation for South monthly averages
acf(monthly_avg_south$avg_eac_cci, main = "Autocorrelation - South Region EAC CCI")

# Autocorrelation for North monthly averages
acf(monthly_avg_north$avg_eac_cci, main = "Autocorrelation - North Region EAC CCI")



df_clean <- df %>% filter(!is.na(South) & !is.na(North))
# Fit linear model: North as a function of South
lm_fit_north <- lm(North ~ South, data = df_clean)

par(mfrow = c(2, 2))  # Arrange plots in a 2x2 grid
plot(lm_fit_north)


# Extract residuals
residuals_north <- resid(lm_fit_north)

# You can also fit South as a function of North if you want both sets of residuals:
lm_fit_south <- lm(South ~ North, data = df_clean)
residuals_south <- resid(lm_fit_south)

# Correlate the residuals
residuals_correlation <- cor(residuals_north, residuals_south, method = "pearson")
residuals_cor_test <- cor.test(residuals_north, residuals_south, method = "pearson")

# Output results
print(paste("Pearson correlation of residuals:", round(residuals_correlation, 3)))
print(residuals_cor_test)



# Create a dataframe for residuals
residuals_df <- data.frame(
  Residuals_North = residuals_north,
  Residuals_South = residuals_south
)

# Plot residuals correlation
ggplot(residuals_df, aes(x = Residuals_South, y = Residuals_North)) +
  geom_point(color = "blue", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  labs(
    title = "Correlation of Residuals (South vs North)",
    x = "Residuals (South ~ North)",
    y = "Residuals (North ~ South)",
    subtitle = paste("Pearson correlation:", round(residuals_correlation, 3),
                     ", p-value = ", signif(residuals_cor_test$p.value, 3))
  ) +
  theme(axis.title = element_text(size = 14),
        plot.title = element_text(size = 18, face = "bold"))

ggsave("output/eac_cci_residuals_correlation.png", width = 8, height = 6, dpi = 300)

# Get fitted values from both models
fitted_north <- fitted(lm_fit_north)   # Predicted North from South
fitted_south <- fitted(lm_fit_south) # Predicted South from North


# Run Pearson correlation
fitted_correlation <- cor(fitted_south, fitted_north, method = "pearson")
fitted_cor_test <- cor.test(fitted_south, fitted_north, method = "pearson")

# Output results
print(paste("Pearson correlation of fitted values:", round(fitted_correlation, 3)))
print(fitted_cor_test)


# Create dataframe for plotting
fitted_df <- data.frame(
  Fitted_North = fitted_north,
  Fitted_South = fitted_south
)

# Plot fitted values against each other
ggplot(fitted_df, aes(x = Fitted_South, y = Fitted_North)) +
  geom_point(color = "blue", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  labs(
    title = "Comparison of Fitted Values: North vs South",
    x = "Fitted South (from North)",
    y = "Fitted North (from South)",
    subtitle = paste("Pearson correlation:", round(fitted_correlation, 3),
                     ", p-value = ", signif(fitted_cor_test$p.value, 3))
  ) +
  theme(axis.title = element_text(size = 14),
        plot.title = element_text(size = 18, face = "bold"))

ggsave("output/eac_cci_fitted_comparison.png", width = 8, height = 6, dpi = 300)


draw(lm_fit)



print(summary(lm_fit))
print(summary(lm_fit_south))




# linear regression of North v South

# Fit a linear model
ns_model <- lm(North ~ South, data = df_clean)

# Summary of the model
summary(ns_model)

# Plot the regression
plot(df_clean$North, df_clean$South, main = "Regression of North EAC CCI on South EAC CCI",
     xlab = "North EAC CCI", ylab = "South EAC CCI", pch = 19)
abline(model, col = "blue", lwd = 2)



plot(ns_model$fitted.values, ns_model$residuals,
     main = "Residuals vs Fitted Values",
     xlab = "Fitted Values",
     ylab = "Residuals",
     pch = 19, col = "darkgreen")
abline(h = 0, col = "red", lwd = 2)

par(mfrow = c(2, 2))  # Arrange plots in a 2x2 grid
plot(ns_model)
