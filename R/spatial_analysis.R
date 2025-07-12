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
