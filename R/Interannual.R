### interannual.R
###
### Analyse interannual difference in the Index
###
### Created: 2025-09-06
### Author: Davina Gifford
### Last updated: 2025-09-06
### Edited by: Davina Gifford

# load libraries

library(tidyverse)

# Get index data

index <- month_data

# Identify seasons for the data

index <- index %>%
  mutate(season = case_when(
    month(trip_month) %in% c(12, 1, 2) ~ "Summer",
    month(trip_month) %in% c(3, 4, 5) ~ "Autumn",
    month(trip_month) %in% c(6, 7, 8) ~ "Winter",
    month(trip_month) %in% c(9, 10, 11) ~ "Spring"
  ))

# compare seasons across years

index_summary <- index %>%
  group_by(year = year(trip_month), season) %>%
  summarise(mean_eac_cci = mean(eac_cci, na.rm = TRUE)) %>%
  ungroup()
print(index_summary)


# Plot seasonal means across years
ggplot(index_summary, aes(x = year, y = mean_eac_cci, color = season)) +
  geom_line() +
  geom_point() +
  labs(title = "Seasonal Mean EAC CCI Across Years",
       x = "Year",
       y = "Mean EAC CCI") 


# add trend line to the seasonal means plot 

ggplot(index_summary, aes(x = year, y = mean_eac_cci, color = season)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_smooth(method = "loess", se = FALSE, linetype = "dashed") +
  labs(title = "Seasonal Mean EAC CCI Across Years with Trend Line",
       x = "Year",
       y = "Mean EAC CCI") 


# facet the plot by season

ggplot(index_summary, aes(x = year, y = mean_eac_cci)) +
  geom_line(color = "blue", linewidth = 1.2) +
  geom_point(size = 2, color = "blue") +
  geom_smooth(method = "loess", se = TRUE, linetype = "dashed", color = "black") +
  labs(title = "Seasonal Mean EAC CCI Across Years by Season",
       x = "Year",
       y = "Mean EAC CCI") +
  facet_wrap(~ season) 



  
# plot seasonal means with error bars

seasonal_summary <- index %>%
  group_by(year = year(trip_month), season) %>%
  summarise(
    mean_eac_cci = mean(eac_cci, na.rm = TRUE),
    sd_eac_cci = sd(eac_cci, na.rm = TRUE),
    n = n(),
    se_eac_cci = sd_eac_cci / sqrt(n)
  ) %>%
  ungroup()

# plot the seasonal means with error bars
ggplot(seasonal_summary, aes(x = year, y = mean_eac_cci, color = season)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean_eac_cci - se_eac_cci, ymax = mean_eac_cci + se_eac_cci), width = 0.2) +
  labs(title = "Seasonal Mean EAC CCI with Standard Error",
       x = "Year",
       y = "Mean EAC CCI") +
  theme_minimal()


# calculate relationships across years

yearly_summary <- index %>%
  group_by(year = year(trip_month)) %>%
  summarise(
    mean_eac_cci = mean(eac_cci, na.rm = TRUE),
    sd_eac_cci = sd(eac_cci, na.rm = TRUE),
    n = n(),
    se_eac_cci = sd_eac_cci / sqrt(n)
  ) %>%
  ungroup()
print(yearly_summary)


# Plot yearly means with error bars
ggplot(yearly_summary, aes(x = year, y = mean_eac_cci)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean_eac_cci - se_eac_cci, ymax = mean_eac_cci + se_eac_cci), width = 0.2) +
  labs(title = "Yearly Mean EAC CCI with Standard Error",
       x = "Year",
       y = "Mean EAC CCI") +
  theme_minimal()


# understand the diferences between seasons
seasonal_diff <- index %>%
  group_by(season) %>%
  summarise(
    mean_eac_cci = mean(eac_cci, na.rm = TRUE),
    sd_eac_cci = sd(eac_cci, na.rm = TRUE),
    n = n(),
    se_eac_cci = sd_eac_cci / sqrt(n)
  ) %>%
  ungroup()
print(seasonal_diff)
# Plot seasonal differences\
ggplot(seasonal_diff, aes(x = season, y = mean_eac_cci, fill = season)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mean_eac_cci - se_eac_cci, ymax = mean_eac_cci + se_eac_cci), width = 0.2) +
  labs(title = "Mean EAC CCI by Season with Standard Error",
       x = "Season",
       y = "Mean EAC CCI") +
  theme_minimal() +
  theme(legend.position = "none")
# ANOVA to test for significant differences between seasons
anova_result <- aov(eac_cci ~ season, data = index)
summary(anova_result)
# Post-hoc test (Tukey's HSD) if ANOVA is significant
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)
# Visualize Tukey's HSD results
plot(tukey_result)
