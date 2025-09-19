### anomaly_comparison.R
###
### Compare the anomalies of the Index with the anomalies from the strength values calculated from the Copernicus Data
###
### File created: 2025-09-16
### Created by: Davina Gifford
### File updated: 2025-09-16
### Updated by: Davina Gifford


# load libraries
library(tidyverse)
library(lubridate)
library(ggplot2)


# Raw index and strength data -----------------------------------------------

index_anomaly <- month_data_t %>% 
  select(trip_month, eac_cci)


strength_anomaly <- strength_with_clim_full_50 %>% 
  select(date, mean_vel)

# Rename columns for clarity
strength_anomaly <- strength_anomaly %>%
  rename(trip_month = date, eac_strength = mean_vel)
# Merge the two datasets on trip_month
combined_anomalies <- index_anomaly %>%
  inner_join(strength_anomaly, by = "trip_month")

# pearson correlation between the two anomalies
corr_anom <- cor(combined_anomalies$eac_cci, combined_anomalies$eac_strength, method = "pearson", use = "complete.obs")
corr_test_anom <- cor.test(combined_anomalies$eac_cci, combined_anomalies$eac_strength, method = "pearson", use = "complete.obs")


# Output results
print(paste("Pearson correlation:", round(corr_anom, 3)))
print(corr_test_anom)


p_val_anom <- corr_test_anom$p.value


# kendall correlation

corr_anom_k <- cor(combined_anomalies$eac_cci, combined_anomalies$eac_strength, method = "kendall", use = "complete.obs")
corr_test_anom_k <- cor.test(combined_anomalies$eac_cci, combined_anomalies$eac_strength, method = "kendall", use = "complete.obs")


# Output results
print(paste("Kendall correlation:", round(corr_anom_k, 3)))
print(corr_test_anom_k)


p_val_anom_k <- corr_test_anom_k$p.value



# plot a scatterplot of the two

ggplot(data = combined_anomalies) +
  geom_point(mapping = aes(x = eac_cci, y = eac_strength)) +
  geom_smooth(mapping = aes(x = eac_cci, y = eac_strength), method = "loess", se = TRUE, color = "blue", linetype = "dashed") +
  labs(
    x = "EAC CCI ",
    y = "Current Strength ",
    title = "Scatterplot of Current Strength vs EAC CCI ",
    subtitle = paste("Pearson correlation:", round(corr_anom, 3), "| p-value:", signif(p_val_anom, 3))
  ) 


# plot the two against trip month

ggplot(data = combined_anomalies, aes(x = trip_month)) +
  geom_line(aes(y = eac_cci, color = "EAC CCI"), size = 1) +
  geom_point(aes(y = eac_cci, color = "EAC CCI"), size = 2) +
  geom_smooth(aes(y = eac_cci, color = "EAC CCI"), method = "loess", se = TRUE, linetype = "dashed") +
  geom_line(aes(y = eac_strength, color = "Current Strength"), size = 1) +
  geom_point(aes(y = eac_strength, color = "Current Strength"), size = 2) +
  geom_smooth(aes(y = eac_strength, color = "Current Strength"), method = "loess", se = TRUE, linetype = "dashed") +
  scale_color_manual(values = c("EAC CCI" = "black", "Current Strength" = "blue")) +
  labs(
    x = "Time",
    y = "Values",
    title = "EAC CCI and Current Strength Over Time",
    color = "Legend"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top"
  )


# plot with the index values



strength_scale <- max(combined_anomalies$eac_cci, na.rm = TRUE) / max(abs(combined_anomalies$eac_strength), na.rm = TRUE)

combined_anomalies <- combined_anomalies %>%
  mutate(str_scale = eac_strength * strength_scale)

combined_anomalies <- combined_anomalies %>%
  mutate(trip_month = as.Date(trip_month))

p <- ggplot(combined_anomalies, aes(x = trip_month)) +
  geom_line(aes(y = eac_cci), colour = "black", linewidth = 1) +
  geom_point(aes(y = eac_cci), size = 3, shape = 21, fill = "black") +
  geom_smooth(aes(y = eac_cci), method = "loess", se = TRUE, color = "black", linetype = "dashed") +
  geom_line(aes(y = eac_strength), colour = "blue", linewidth = 1) +
  geom_point(aes(y = eac_strength), size = 3, shape = 21, fill = "blue") +
  geom_smooth(aes(y = eac_strength), method = "loess", se = TRUE, color = "blue", linetype = "dashed") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(name = "EAC copepod composition index",
                     sec.axis = sec_axis(~ .,  name = "Monthy mean strength")) +
  theme(
    axis.title.y = element_text(size = 14, color = "black"),
    axis.title.y.right = element_text(size = 14, color = "blue"),
    plot.title = element_text(size = 18, face = "bold")
  ) +
  labs(
    x = "Time",
    title = "EAC copepod composition index & Monthly mean strength ")



ggsave(file.path("output", "eac_cci_str_combined.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)


# anomalies --------------------

index_anomaly_2 <- month_data_t %>% 
  select(trip_month, eac_cci_clim)


strength_anomaly2 <- strength_with_clim_full_50 %>% 
  select(date, vel_anomaly)

# Rename columns for clarity
strength_anomaly2 <- strength_anomaly2 %>%
  rename(trip_month = date, strength_anom = vel_anomaly)

# Merge the two datasets on trip_month
combined_anomalies2 <- index_anomaly_2 %>%
  full_join(strength_anomaly2, by = "trip_month")

# pearson correlation between the two anomalies
corr_anom2 <- cor(combined_anomalies2$eac_cci_clim, combined_anomalies2$strength_anom, method = "pearson", use = "complete.obs")
corr_test_anom2 <- cor.test(combined_anomalies2$eac_cci_clim, combined_anomalies2$strength_anom, method = "pearson", use = "complete.obs")


# Output results
print(paste("Pearson correlation:", round(corr_anom2, 3)))
print(corr_test_anom2)


p_val_anom2 <- corr_test_anom2$p.value


# kendall correlation between the two anomalies
corr_anom2_k <- cor(combined_anomalies2$eac_cci_clim, combined_anomalies2$strength_anom, method = "kendall", use = "complete.obs")
corr_test_anom2_k <- cor.test(combined_anomalies2$eac_cci_clim, combined_anomalies2$strength_anom, method = "kendall", use = "complete.obs")


# Output results
print(paste("Kendall correlation:", round(corr_anom2_k, 3)))
print(corr_test_anom2_k)


p_val_anom2_k <- corr_test_anom2_k$p.value


# plot the two 



strength_anom_scale <- max(combined_anomalies2$eac_cci_clim, na.rm = TRUE) /
  max(abs(combined_anomalies2$strength_anom), na.rm = TRUE)


combined_anomalies2 <- combined_anomalies2 %>%
  mutate(str_scale = strength_anom * strength_anom_scale,
         trip_month = as.Date(trip_month)) %>%
  filter(
    !is.na(eac_cci_clim),
    !is.na(str_scale),
    is.finite(eac_cci_clim),
    is.finite(str_scale)
  )


p <- ggplot(combined_anomalies2, aes(x = trip_month)) +
  geom_line(aes(y = eac_cci_clim), colour = "black", linewidth = 1) +
  geom_point(aes(y = eac_cci_clim), size = 3, shape = 21, fill = "black") +
  geom_smooth(aes(y = eac_cci_clim), method = "loess", se = TRUE, color = "black", linetype = "dashed") +
  geom_line(aes(y = str_scale), colour = "blue", linewidth = 1) +
  geom_point(aes(y = str_scale), size = 3, shape = 21, fill = "blue") +
  geom_smooth(aes(y = str_scale), method = "loess", se = TRUE, color = "blue", linetype = "dashed") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(name = "EAC copepod composition index",
                     sec.axis = sec_axis(~ ./ strength_anom_scale,  name = "Monthly mean strength")) +
  theme(
    axis.title.y = element_text(size = 14, color = "black"),
    axis.title.y.right = element_text(size = 14, color = "blue"),
    plot.title = element_text(size = 18, face = "bold")
  ) +
  labs(
    x = "Time",
    title = "EAC copepod composition index & Monthly mean strength ")



ggsave(file.path("output", "eac_cci_str_combined.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)

# try to get the plotting fixed


# Load required package
library(zoo)

# Ensure data is sorted by date
combined_anomalies2 <- combined_anomalies2 %>%
  arrange(trip_month)

# Interpolate missing values for eac_cci_clim and strength_anom
combined_anomalies2 <- combined_anomalies2 %>%
  mutate(
    eac_cci_clim_interp = na.approx(eac_cci_clim, x = trip_month, na.rm = FALSE),
    strength_anom_interp = na.approx(strength_anom, x = trip_month, na.rm = FALSE)
  )

strength_anom_scale <- max(combined_anomalies2$eac_cci_clim_interp, na.rm = TRUE) /
  max(abs(combined_anomalies2$strength_anom_interp), na.rm = TRUE)

combined_anomalies2 <- combined_anomalies2 %>%
  mutate(
    str_scale = strength_anom_interp * strength_anom_scale,
    trip_month = as.Date(trip_month)
  )

p <- ggplot(combined_anomalies2, aes(x = trip_month)) +
  geom_line(aes(y = eac_cci_clim_interp), colour = "black", linewidth = 1, na.rm = TRUE) +
  geom_point(aes(y = eac_cci_clim_interp), size = 3, shape = 21, fill = "black", na.rm = TRUE) +
  geom_smooth(aes(y = eac_cci_clim_interp), method = "loess", se = TRUE, color = "black", linetype = "dashed", na.rm = TRUE) +
  geom_line(aes(y = str_scale), colour = "blue", linewidth = 1, na.rm = TRUE) +
  geom_point(aes(y = str_scale), size = 3, shape = 21, fill = "blue", na.rm = TRUE) +
  geom_smooth(aes(y = str_scale), method = "loess", se = TRUE, color = "blue", linetype = "dashed", na.rm = TRUE) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(
    name = "EAC copepod composition index",
    sec.axis = sec_axis(~ . / strength_anom_scale, name = "Monthly mean strength")
  ) +
  theme(
    axis.title.y = element_text(size = 14, color = "black"),
    axis.title.y.right = element_text(size = 14, color = "blue"),
    plot.title = element_text(size = 18, face = "bold")
  ) +
  labs(
    x = "Time",
    title = "EAC copepod composition index & Monthly mean strength"
  )

ggsave(file.path("output", "eac_cci_str_combined_interp.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = "png")



# do anova 

anova_data <- combined_anomalies2 %>%
  select(trip_month, eac_cci_clim, strength_anom) %>%
  filter(!is.na(eac_cci_clim) & !is.na(strength_anom))
anova_data <- anova_data %>%
  mutate(month = month(trip_month, label = TRUE))
anova_data$month <- as.factor(anova_data$month)
anova_data$year <- year(anova_data$trip_month)
anova_data$year <- as.factor(anova_data$year)
anova_data <- anova_data %>%
  filter(!is.na(month) & !is.na(year))
anova_data <- na.omit(anova_data)

anova_result <- aov(eac_cci_clim ~ strength_anom + month + year, data = anova_data)
summary(anova_result)

# Residuals vs Fitted
plot(anova_result$fitted.values, anova_result$residuals,
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted Values",
     pch = 20, col = "blue")
abline(h = 0, col = "red", lty = 2)

# QQ Plot for normality of residuals
qqnorm(anova_result$residuals, main = "QQ Plot of Residuals")
qqline(anova_result$residuals, col = "red")


# save to rds

saveRDS(combined_anomalies2, file = file.path("var", "combined_anomalies.rds"))
