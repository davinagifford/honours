### anomaly_analysis_south.R
###
### Exploring the anomalies and how they differ for the samples south of the separation zone. 
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



anomaly_south <- tibble(
  trip_month = month_data_t$trip_month,
  observed_eac_cci = month_data_t$eac_cci,
  climatology = month_data_t$eac_cci_clim,
  anomaly = month_data_t$eac_cci - month_data_t$eac_cci_clim
)


anomaly_south$trip_month <- as.Date(anomaly_south$trip_month)


# Create the plot with a trend line
ggplot(anomaly_south, aes(x = trip_month, y = anomaly)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "darkred", size = 2) +
  geom_smooth(method = "loess", se = TRUE, color = "black", linetype = "dashed") +
  labs(
    title = "EAC CCI Anomaly Over Time with Trend Line - south",
    x = "Trip Month",
    y = "EAC CCI Anomaly"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# save image 
# ggsave("EAC_CCI_Anomaly_Over_Time.png", width = 1200 / 96, height = 600 / 96, dpi = 96)


# Set threshold for strong anomaly (z-score)
threshold <- 1 # Adjust as needed

# Calculate anomaly strength and classify direction
anomaly_south <- anomaly_south %>%
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
ggplot(anomaly_south, aes(x = trip_month, y = anomaly_z)) +
  geom_line() +
  geom_point(aes(color = is_strong), size = 2) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = paste("EAC Climate Anomalies (z-score) with Strong Events Highlighted (Threshold =", threshold, ") - south"),
    y = "Anomaly (z-score)", x = "Month", color = "Strong Anomaly"
  ) +
  theme_minimal() 
#+
#  facet_wrap(~ season)

# Classify direction of strong anomalies
strong_anomalies_n <- anomaly_south %>%
  filter(is_strong)

# Count how many are positive vs negative
table(strong_anomalies_n$direction)

# Visualise direction of strong anomalies
ggplot(strong_anomalies_n, aes(x = direction)) +
  geom_bar(fill = "steelblue") +
  labs(
    title = "Direction of Strong Anomalies - south",
    x = "Direction", y = "Count"
  ) +
  theme_minimal()

# Run ANOVA
anova_model_n <- aov(anomaly_strength ~ season, data = anomaly_south)
summary(anova_model_n)

# Full LOESS-smoothed anomaly strength plot
ggplot(anomaly_south, aes(x = trip_month, y = anomaly_strength)) +
  geom_point(color = "black", alpha = 0.6) +
  geom_smooth(method = "loess", span = 0.2, se = FALSE, color = "green", linewidth = 1.2) +
  labs(
    title = "LOESS Smoothing of Anomaly Strength Over Time (z-score) - south",
    x = "Month",
    y = "Anomaly Strength (|z-score|)"
  ) +
  theme_minimal()

# Plot with separate LOESS lines by direction
ggplot(anomaly_south, aes(x = trip_month, y = anomaly_strength, color = direction)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", span = 0.3, se = FALSE, linewidth = 1, linetype = "dashed") +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  labs(
    title = "LOESS Smoothing of Anomaly Strength by Direction (z-score) - south",
    x = "Year",
    y = "Anomaly Strength (|z-score|)"
  ) +
  scale_color_manual(values = c("Positive" = "red", "Negative" = "blue")) +
  theme_minimal() 
#+
#  facet_wrap(~ season)


anomaly_south$season <- factor(anomaly_south$season, levels = c("Winter", "Spring", "Summer", "Autumn"))

# seasonal breakdown of anomalies
# Basic boxplot
ggplot(anomaly_south, aes(x = season, y = anomaly, fill = season)) +
  geom_boxplot(alpha = 0.6) +
  labs(
    title = "EAC CCI Anomalies by Season - south",
    x = "Season",
    y = "EAC CCI Anomaly"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )
