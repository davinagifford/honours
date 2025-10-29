### spatial_analysis.R
###
### Analysis of the differences between north and south of the EAC.
###
### Created: 2025-07-04
### Author: Davina Gifford
### Last updated: 2025-10-30
### Edited by: Davina Gifford

# load libraries
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyverse)

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
correlation_ns <- cor(south_values, north_values, method = "pearson", use = "complete.obs")
cor_test_ns <- cor.test(south_values, north_values, method = "pearson", use = "complete.obs")

# Output results
print(paste("Pearson correlation:", round(correlation_ns, 3)))
print(cor_test_ns)


p_val_ns <- cor_test_ns$p.value


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
    subtitle = paste("Pearson correlation:", round(correlation_ns, 3),
                     ", p-value = ", signif(p_val_ns, 3))
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









# Residuals  --------------------------------------------------------------
# linear regression of North v South

# Fit a linear model
ns_model <- lm(North ~ South, data = df_clean)

# Summary of the model
summary(ns_model)



residuals_ns <- resid(ns_model)
predicted_ns <- fitted(ns_model)



plot(df_clean$South, ns_model$residuals,
     main = "Residuals vs Fitted Values",
     xlab = "Fitted Values",
     ylab = "Residuals",
     pch = 19, col = "darkgreen")
abline(h = 0, col = "red", lwd = 2)



# scatterplot of north v south -----------

# Pivot longer for easier plotting


df_long <- df %>%
  pivot_longer(cols = c(South, North),
               names_to = "Region",
               values_to = "value")

df_long$Month <- factor(df_long$Month, levels = month.abb)  # if Month is abbreviated

df_long <- df_long %>%
  arrange(Region, Month)


df_long$Month_num <- match(df_long$Month, month.abb)



p <- ggplot(df_long, aes(x = Month_num, y = value, color = Region)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3, shape = 21, fill = "white") +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_color_manual(values = c("North" = "blue", "South" = "red"),
                     labels = c("North", "South")) +
  labs(
    x = "Month",
    y = "Region value",
    color = "Region",
    title = "Comparison of North and South EAC CCI values",
    subtitle = paste("Pearson correlation:", round(correlation_ns, 3), "| p-value:", signif(p_val_ns, 3))
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 12),
  )

ggsave(file.path("output", "north-south-compare.png"),
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




# correlate the raw values not monthly averages --------


north_south_data <- month_data_north %>% 
  rename(north_cci = eac_cci) %>% 
  left_join(month_data_south, by = "trip_month") %>% 
  rename(south_cci = eac_cci)

north_south_data$trip_month <- as.Date(north_south_data$trip_month)
  


correlation_raw <- cor(north_south_data$north_cci, north_south_data$south_cci, method = "pearson", use = "complete.obs")
cor_test_raw <- cor.test(north_south_data$north_cci, north_south_data$south_cci, method = "pearson", use = "complete.obs")

# Output results
print(paste("Pearson correlation (raw):", round(correlation_raw, 3)))
print(cor_test_raw)

# plot the two 





p <- ggplot(north_south_data, aes(x = trip_month)) +
  geom_line(aes(y = north_cci), colour = "black", linewidth = 1) +
  geom_point(aes(y = north_cci), size = 3, shape = 21, fill = "black") +
  geom_smooth(aes(y = north_cci), method = "loess", se = TRUE, color = "black", linetype = "dashed") +
  geom_line(aes(y = south_cci), colour = "blue", linewidth = 1) +
  geom_point(aes(y = south_cci), size = 3, shape = 21, fill = "blue") +
  geom_smooth(aes(y = south_cci), method = "loess", se = TRUE, color = "blue", linetype = "dashed") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(name = "EAC copepod composition index (North)",
                     sec.axis = sec_axis(~ .,  name = "EAC copepod composition index (South)")) +
  theme(
    axis.title.y = element_text(size = 14, color = "black"),
    axis.title.y.right = element_text(size = 14, color = "blue"),
    plot.title = element_text(size = 18, face = "bold")
  ) +
  labs(
    x = "Time",
    title = "North and South EAC copepod composition index  ")



ggsave(file.path("output", "eac_cci_NS_combined.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)

# do anova  -----

anova_result_ns <- aov(north_cci ~ south_cci, data = north_south_data)

summary(anova_result_ns)

# Residuals vs Fitted
plot(anova_result_ns$fitted.values, anova_result_ns$residuals,
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted Values",
     pch = 20, col = "blue")
abline(h = 0, col = "red", lty = 2)
# Q-Q Plot

qqnorm(anova_result_ns$residuals, main = "Q-Q Plot of Residuals")
qqline(anova_result_ns$residuals, col = "red", lty = 2)



# compare north south anomalies -----------------

north_anomaly <- month_data_t_n %>% 
  select(trip_month, eac_cci_clim) %>% 
  rename(north_anom = eac_cci_clim)

north_anomaly <- north_anomaly %>% 
  mutate(north_anom = north_anom - month_data_t_n$eac_cci)

south_anomaly <- month_data_t_s %>% 
  select(trip_month, eac_cci_clim) %>% 
  rename(south_anom = eac_cci_clim)

south_anomaly <- south_anomaly %>% 
  mutate(south_anom = south_anom - month_data_t_s$eac_cci)

ns_anom_combined <- north_anomaly %>% 
  left_join(south_anomaly, by = "trip_month")

ns_anom_combined$trip_month <- as.Date(ns_anom_combined$trip_month)


# pearson correlation between the two anomalies
corr_anom_ns <- cor(ns_anom_combined$north_anom, ns_anom_combined$south_anom, method = "pearson", use = "complete.obs")
corr_test_anom_ns <- cor.test(ns_anom_combined$north_anom, ns_anom_combined$south_anom, method = "pearson", use = "complete.obs")


# Output results
print(paste("Pearson correlation:", round(corr_anom_ns, 3)))
print(corr_test_anom_ns)


p_val_anom_ns <- corr_test_anom_ns$p.value

# kendall correlation between the two anomalies
corr_anom_ns_k <- cor(ns_anom_combined$north_anom, ns_anom_combined$south_anom, method = "kendall", use = "complete.obs")
corr_test_anom_ns_k <- cor.test(ns_anom_combined$north_anom, ns_anom_combined$south_anom, method = "kendall", use = "complete.obs")


# Output results
print(paste("kendall correlation:", round(corr_anom_ns_k, 3)))
print(corr_test_anom_ns_k)


p_val_anom_ns_k <- corr_test_anom_ns_k$p.value



# plot the two anomalies

p <- ggplot(ns_anom_combined, aes(x = trip_month)) +
  geom_line(aes(y = north_anom), colour = "black", linewidth = 1) +
  geom_point(aes(y = north_anom), size = 3, shape = 21, fill = "black") +
  geom_smooth(aes(y = north_anom), method = "loess", se = TRUE, color = "black", linetype = "dashed") +
  geom_line(aes(y = south_anom), colour = "blue", linewidth = 1) +
  geom_point(aes(y = south_anom), size = 3, shape = 21, fill = "blue") +
  geom_smooth(aes(y = south_anom), method = "loess", se = TRUE, color = "blue", linetype = "dashed") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(name = "EAC copepod composition index - North",
                     sec.axis = sec_axis(~ .,  name = "EAC copepod composition index - South")) +
  theme(
    axis.title.y = element_text(size = 14, color = "black"),
    axis.title.y.right = element_text(size = 14, color = "blue"),
    plot.title = element_text(size = 18, face = "bold")
  ) +
  labs(
    x = "Time",
    title = "North and South Anomalies of EAC copepod composition index ")



ggsave(file.path("output", "eac_cci_ns_anom_combined.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)
# scatterplot of the two anomalies
p <- ggplot(data = ns_anom_combined) +
  geom_point(mapping = aes(x = north_anom, y = south_anom)) +
  geom_smooth(mapping = aes(x = north_anom, y = south_anom), method = "lm", se = TRUE, color = "blue", linetype = "dashed") +
  labs(
    x = "North Anomaly",
    y = "South Anomaly",
    title = "Scatterplot of North vs South EAC CCI Anomalies",
    subtitle = paste("Pearson correlation:", round(corr_anom_ns, 3), "| p-value:", signif(p_val_anom_ns, 3))
  )
ggsave(file.path("output", "eac_cci_ns_anom_scatter.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)


# plot residuals

ns_anom_model <- lm(south_anom ~ north_anom, data = ns_anom_combined)

summary(ns_anom_model)

residuals_ns_anom <- resid(ns_anom_model)
predicted_ns_anom <- fitted(ns_anom_model)

plot(ns_anom_combined$north_anom, ns_anom_model$residuals,
     main = "Residuals vs Fitted Values",
     xlab = "Fitted Values",
     ylab = "Residuals",
     pch = 19, col = "darkgreen")
abline(h = 0, col = "red", lwd = 2)


# anova of anomaly values
anova_result_ns_anom <- aov(north_anom ~ south_anom, data = ns_anom_combined)

summary(anova_result_ns_anom)

# Residuals vs Fitted
plot(anova_result_ns_anom$fitted.values, anova_result_ns_anom$residuals,
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted Values",
     pch = 20, col = "blue")
abline(h = 0, col = "red", lty = 2)
# Q-Q Plot

qqnorm(anova_result_ns_anom$residuals, main = "Q-Q Plot of Residuals")
qqline(anova_result_ns_anom$residuals, col = "red", lty = 2)
