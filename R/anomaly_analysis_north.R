### anomaly_analysis.R
###
### Exploring the anomalies and how they differ
###
### Created: 2025-08-08
### Author: Davina Gifford
### Last updated: 2025-08-08
### Edited by: Davina Gifford


# Load required libraries
library(tidyverse)
library(ggplot2)
library(readxl)
library(scales)
library(dplyr)

# EAC CCI Anomaly over time -----------------------------------------------



# get anomaly data

anomaly_north <- tibble(
  trip_month = month_data_t$trip_month,
  eac_cci_clim = month_data_t$eac_cci_clim)


month_data_t$trip_month <- as.Date(month_data_t$trip_month)


# Create the plot with a trend line
ggplot(month_data_t, aes(x = trip_month, y = eac_cci_clim)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "darkred", size = 2) +
  geom_smooth(method = "loess", se = TRUE, color = "black", linetype = "dashed") +
  labs(
    title = "EAC CCI Clim Over Time with Trend Line - North",
    x = "Trip Month",
    y = "EAC CCI Clim"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Convert to time series object (assuming monthly data)
ts_data <- ts(month_data_t$eac_cci_clim, start = c(2010, 1), frequency = 12)

# Decompose the time series
decomposed <- decompose(ts_data)

# Plot the decomposition
plot(decomposed)


ts_data <- ts(month_data_t$eac_cci_clim, start = c(2010, 1), frequency = 12)
decomposed_stl <- stl(ts_data, s.window = "periodic")
plot(decomposed_stl)






# Load your data
data <- read_excel("data/Anomaly_north.xlsx")

# Convert date column
data$trip_month <- as.Date(data$trip_month)

# Create a gradient-colored line plot
ggplot(data, aes(x = trip_month, y = eac_cci_clim, color = eac_cci_clim)) +
  geom_point(size = 3) +
  geom_line() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(title = "Gradient Plot of eac_cci_clim",
       x = "Date",
       y = "Value",
       color = "Gradient")



# Create anomaly strength
df <- data %>%
  mutate(
    anomaly_strength = abs(eac_cci_clim),
    is_strong = anomaly_strength >= 0.07
  )

# Plot
ggplot(df, aes(x = trip_month, y = eac_cci_clim)) +
  geom_line() +
  geom_point(aes(color = is_strong), size = 2) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "EAC Climate Anomalies with Strong Events Highlighted - North",
       y = "Anomaly", x = "Month", color = "Strong Anomaly") +
  theme_minimal()


# Classify direction of strong anomalies
strong_anomalies <- df %>%
  filter(anomaly_strength >= 0.07) %>%
  mutate(direction = ifelse(eac_cci_clim > 0, "Positive", "Negative"))

# Count how many are positive vs negative
table(strong_anomalies$direction)

# Optional: Visualise
ggplot(strong_anomalies, aes(x = direction)) +
  geom_bar(fill = "steelblue") +
  labs(title = "Direction of Strong Anomalies - North",
       x = "Direction", y = "Count") +
  theme_minimal()


# Run ANOVA
anova_model <- aov(anomaly_strength ~ Season, data = df)
summary(anova_model)





# Ensure date format and compute anomaly strength
df <- df %>%
  mutate(
    trip_month = as.Date(trip_month),
    anomaly_strength = abs(eac_cci_clim)
  )

# Full LOESS-smoothed anomaly strength plot
ggplot(df, aes(x = trip_month, y = anomaly_strength)) +
  geom_point(color = "black", alpha = 0.6) +
  geom_smooth(method = "loess", span = 0.2, se = FALSE, color = "green", linewidth = 1.2) +
  labs(
    title = "LOESS Smoothing of Anomaly Strength Over Time - North",
    x = "Month",
    y = "Anomaly Strength (|Value|)"
  ) +
  theme_minimal()

# Classify direction
df <- df %>%
  mutate(direction = ifelse(eac_cci_clim > 0, "Positive", "Negative"))


# Extract year data
df <- df %>% 
  separate(trip_month, c("year", "month", "day"), sep = "-", remove = FALSE)

# Plot with separate LOESS lines by direction
ggplot(df, aes(x = trip_month, y = anomaly_strength, color = direction)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", span = 0.3, se = FALSE, linewidth = 1, linetype = "dashed") +
  scale_x_date(
    date_breaks = "1 year",       # Show every year
    date_labels = "%Y"            # Format as 4-digit year
  ) +
  labs(
    title = "LOESS Smoothing of Anomaly Strength by Direction - North",
    x = "Year",
    y = "Anomaly Strength (|Value|)"
  ) +
  scale_color_manual(values = c("Positive" = "red", "Negative" = "blue")) +
  theme_minimal()



