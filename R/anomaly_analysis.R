### anomaly_analysis.R
###
### Exploring the anomalies and how they differ
###
### Created: 2025-08-08
### Author: Davina Gifford
### Last updated: 2025-08-13
### Edited by: Davina Gifford


# Load required libraries
library(tidyverse)
library(ggplot2)
library(readxl)
library(scales)
library(dplyr)

# EAC CCI Anomaly over time -----------------------------------------------



anomaly <- tibble(
  trip_month = month_data_t$trip_month,
  observed_eac_cci = month_data_t$eac_cci,
  climatology = month_data_t$eac_cci_clim,
  anomaly = month_data_t$eac_cci - month_data_t$eac_cci_clim
)


anomaly$trip_month <- as.Date(anomaly$trip_month)


# Create the plot with a trend line
p <- ggplot(anomaly, aes(x = trip_month, y = anomaly)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "darkred", size = 2) +
  geom_smooth(method = "loess", se = TRUE, color = "black", linetype = "dashed") +
  labs(
    title = "EAC CCI Anomaly Over Time with Trend Line",
    x = "Trip Month",
    y = "EAC CCI Anomaly"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# save image 
ggsave(file.path("output", "cci-anomaly-full.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)

# Set threshold for strong anomaly (z-score)
threshold <- 1 # Adjust as needed

# Calculate anomaly strength and classify direction
anomaly <- anomaly %>%
  mutate(
    anomaly_z = (anomaly - mean(anomaly, na.rm = TRUE)) / sd(anomaly, na.rm = TRUE),
    anomaly_strength = abs(anomaly_z),
    is_strong = anomaly_strength >= threshold,
    direction = ifelse(anomaly > 0, "Positive", "Negative")
  ) %>%
  separate(trip_month, c("year", "month", "day"), sep = "-", remove = FALSE) %>% 
  mutate(season = case_when(
    month %in% c("12", "01", "02") ~ "Summer",
    month %in% c("03", "04", "05") ~ "Autumn",
    month %in% c("06", "07", "08") ~ "Winter",
    month %in% c("09", "10", "11") ~ "Spring"
  ))


# Plot: EAC Climate Anomalies with Strong Events Highlighted
ggplot(anomaly, aes(x = trip_month, y = anomaly_z)) +
  geom_line() +
  geom_point(aes(color = is_strong), size = 2) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = paste("EAC Climate Anomalies (z-score) with Strong Events Highlighted (Threshold =", threshold, ")"),
    y = "Anomaly (z-score)", x = "Month", color = "Strong Anomaly"
  ) +
  theme_minimal() +
  facet_wrap(~ season)

# Classify direction of strong anomalies
strong_anomalies <- anomaly %>%
  filter(is_strong)

# Count how many are positive vs negative
table(strong_anomalies$direction)

# Visualise direction of strong anomalies
ggplot(strong_anomalies, aes(x = direction)) +
  geom_bar(fill = "steelblue") +
  labs(
    title = "Direction of Strong Anomalies",
    x = "Direction", y = "Count"
  ) +
  theme_minimal()

# Run ANOVA
anova_model <- aov(anomaly_strength ~ season, data = anomaly)
summary(anova_model)

# Full LOESS-smoothed anomaly strength plot
ggplot(anomaly, aes(x = trip_month, y = anomaly_strength)) +
  geom_point(color = "black", alpha = 0.6) +
  geom_smooth(method = "loess", span = 0.2, se = FALSE, color = "green", linewidth = 1.2) +
  labs(
    title = "LOESS Smoothing of Anomaly Strength Over Time (z-score)",
    x = "Month",
    y = "Anomaly Strength (|z-score|)"
  ) +
  theme_minimal()

# Plot with separate LOESS lines by direction
ggplot(anomaly, aes(x = trip_month, y = anomaly_strength, color = direction)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", span = 0.3, se = FALSE, linewidth = 1, linetype = "dashed") +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  labs(
    title = "LOESS Smoothing of Anomaly Strength by Direction (z-score)",
    x = "Year",
    y = "Anomaly Strength (|z-score|)"
  ) +
  scale_color_manual(values = c("Positive" = "red", "Negative" = "blue")) +
  theme_minimal() +
  facet_wrap(~ season)


anomaly$season <- factor(anomaly$season, levels = c("Winter", "Spring", "Summer", "Autumn"))

# seasonal breakdown of anomalies
# Basic boxplot
ggplot(anomaly, aes(x = season, y = anomaly, fill = season)) +
  geom_boxplot(alpha = 0.6) +
  labs(
    title = "EAC CCI Anomalies by Season",
    x = "Season",
    y = "EAC CCI Anomaly"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )





# 1. Calculate temperature anomalies (remove climatology)
# Join temp_tbl with its climatology by day-of-year
temp_anom <- temp_tbl %>%
  mutate(doy = yday(date)) %>%
  left_join(
    climatology_str %>% select(doy, temp_clim = str_clim), # replace with your temp climatology if you have one
    by = "doy"
  ) %>%
  mutate(temp_anomaly = mean_temp - temp_clim)

# 2. Calculate EAC index anomalies (remove climatology)
# Assuming you have a climatology for the EAC index (e.g., climatology$eac_cci)
eac_anom <- anomaly %>%
  select(date = trip_month, observed_eac_cci) %>%
  mutate(doy = yday(date)) %>%
  left_join(
    climatology %>% select(sample_time, eac_cci_clim = eac_cci) %>% mutate(doy = yday(sample_time)),
    by = "doy"
  ) %>%
  mutate(eac_cci_anomaly = observed_eac_cci - eac_cci_clim)

# 3. Join the anomalies by date
anom_compare <- temp_anom %>%
  select(date, temp_anomaly) %>%
  left_join(eac_anom %>% select(date, eac_cci_anomaly), by = "date") %>%
  drop_na()

# 4. Correlation and plot
cor_result <- cor.test(anom_compare$temp_anomaly, anom_compare$eac_cci_anomaly, method = "pearson")
print(cor_result)

# Plot
ggplot(anom_compare, aes(x = temp_anomaly, y = eac_cci_anomaly)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    x = "Temperature Anomaly (°C)",
    y = "EAC Index Anomaly",
    title = "Residual Relationship: Temperature vs EAC Index"
  ) +
  annotate("text", 
           x = min(anom_compare$temp_anomaly, na.rm = TRUE), 
           y = max(anom_compare$eac_cci_anomaly, na.rm = TRUE), 
           label = sprintf("r = %.2f\np = %.3g", cor_result$estimate, cor_result$p.value),
           hjust = 0, vjust = 1, size = 5)

# ...existing code...