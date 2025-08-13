### spatial_analysis.R
###
### Analysis of the differences between north and south of the EAC.
###
### Created: 2025-07-04
### Author: Davina Gifford
### Last updated: 2025-07-04
### Edited by: Davina Gifford

# load libraries
library(dplyr)
library(lubridate)
library(ggplot2)

#Calculate monthly averages for south
monthly_avg_south <- south_data %>%
  mutate(month = month(trip_month, label = TRUE)) %>%  # Extract month as a factor with labels
  group_by(month) %>%
  summarise(avg_eac_cci = mean(eac_cci, na.rm = TRUE), .groups = "drop")

print(monthly_avg_south)

#Calculate monthly averages for south
monthly_avg_north <- north_data %>%
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


df_clean <- df %>% filter(!is.na(South) & !is.na(North))
# Fit linear model: North as a function of South
lm_fit <- lm(North ~ South, data = df)

# Extract residuals
residuals_north <- resid(lm_fit)

# You can also fit South as a function of North if you want both sets of residuals:
lm_fit_south <- lm(South ~ North, data = df)
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
