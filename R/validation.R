### validation.R
###
### compare output of index with the original index
###
### Created: 2025-08-19
### Author: Davina Gifford
### Last updated: 2025-08-19
### Edited by: Davina Gifford


library(tidyverse)
library(ncdf4)
library(tidync)

# read mooring data


# in tidy 

datafile <- "D:/HONOURS/DavinaG_2025_Honours/data/EAC_filled-daily-distance-depth-gridded-product_20120401-20220727.nc"
if (!tolower(Sys.info()[["sysname"]]) == "sunos") {
  tnc <- tidync(datafile)
  print(tnc)
}

# explore variables
variables <- hyper_vars(tnc)
print(variables)

# extract velocity data
data_tbl <- tnc %>%
  hyper_tibble(select_var = "VCUR")
print(data_tbl)



raw_nc_data <- tnc %>%
  hyper_tibble(select_var = "VCUR")

# look at the data before doing any mutates
raw_nc_data <- raw_nc_data %>% 
  filter(DEPTH < 1000) %>%
  filter(LONGITUDE > 153.9) %>% 
  filter(LONGITUDE < 155.0)

# filter data to desired depth and longitude range. Create monthly average velocity data
data_tbl <- data_tbl %>%
  filter(DEPTH < 1000) %>%
  filter(LONGITUDE > 153.9) %>% 
  filter(LONGITUDE < 155.0) %>% 
  mutate(
    date = as.Date(sub("T.*", "", TIME)),
    year = year(date),
    month = month(date)
  ) %>%
  group_by(year, month) %>%
  summarise(mean_vcur = mean(VCUR, na.rm = TRUE), .groups = "drop") %>% 
  mutate(date = as.Date(paste(year, month, "01", sep = "-"))) 

# convert velocity to positive strength value
data_tbl <- data_tbl %>%
  mutate(mean_vcur = abs(mean_vcur))


# plot the data
p <- ggplot(data_tbl, aes(x = date, y = mean_vcur)) +
  geom_line(linewidth = 1, colour = "black") +
  geom_point(size = 3, shape = 21, fill = "black") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  geom_smooth(method = "loess", se = TRUE, color = "black", linetype = "dashed") +
  labs(
    x = "Time",
    y = expression("Monthly Mean strength (m/s"^{-1}*")"),
    title = "Mean Strength Over Time (100m depth)"
  ) +
  theme(axis.title = element_text(size = 14), 
        plot.title = element_text(size = 18, face = "bold"))

ggsave(file.path("output", "mean-velocity_100m.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)


# extract temperature data
temp_tbl <- tnc %>%
  hyper_tibble(select_var = "TEMP")
print(temp_tbl)

# filter data to desired depth and longitude range. Create monthly average velocity data
temp_tbl <- temp_tbl %>%
  filter(DEPTH < 1000) %>%
  filter(LONGITUDE > 153.9) %>% 
  filter(LONGITUDE < 155.0) %>% 
  mutate(
    date = as.Date(sub("T.*", "", TIME)),
    year = year(date),
    month = month(date)
  ) %>%
  group_by(year, month) %>%
  summarise(mean_temp = mean(TEMP, na.rm = TRUE), .groups = "drop") %>% 
  mutate(date = as.Date(paste(year, month, "01", sep = "-"))) 



# plot the data
p <- ggplot(temp_tbl, aes(x = date, y = mean_temp)) +
  geom_line(linewidth = 1, colour = "black") +
  geom_point(size = 3, shape = 21, fill = "black") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  geom_smooth(method = "loess", se = TRUE, color = "black", linetype = "dashed") +
  labs(
    x = "Time",
    y = expression("Monthly Mean Temperature"),
    title = "Mean Temperature Over Time (100m depth)"
  ) +
  theme(axis.title = element_text(size = 14), 
        plot.title = element_text(size = 18, face = "bold"))

ggsave(file.path("output", "mean-temp_100m.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)


# is there a relationship between velocity and temperature? 

vel_temp <- data_tbl %>%
  left_join(temp_tbl, by = "date") %>%
  select(date, mean_vcur, mean_temp)

# Rescale for plotting (adjust factor as needed)
temp_scale_factor <- max(vel_temp$mean_vcur, na.rm = TRUE) / max(abs(vel_temp$mean_temp), na.rm = TRUE)
vel_temp <- vel_temp %>%
  mutate(temp_scaled = mean_temp * temp_scale_factor)

# Plot with dual y-axis
p <- ggplot(vel_temp, aes(x = date)) +
  geom_line(aes(y = mean_vcur), colour = "black", linewidth = 1) +
  geom_point(aes(y = mean_vcur), size = 3, shape = 21, fill = "black") +
  geom_smooth(aes(y = mean_vcur), method = "loess", se = TRUE, color = "black", linetype = "dashed") +
  geom_line(aes(y = temp_scaled), colour = "blue", linewidth = 1) +
  geom_point(aes(y = temp_scaled), size = 3, shape = 21, fill = "blue") +
  geom_smooth(aes(y = temp_scaled), method = "loess", se = TRUE, color = "blue", linetype = "dashed") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(
    name = "Mean Monthly strength",
    sec.axis = sec_axis(~ . / temp_scale_factor, name = "Mean monthly Temperature")
  ) +
  theme(
    axis.title.y = element_text(size = 14, color = "black"),
    axis.title.y.right = element_text(size = 14, color = "blue"),
    plot.title = element_text(size = 18, face = "bold")
  ) +
  labs(
    x = "Time",
    title = "Mean monthly temperature and mean monthly strength (100m depth)"
  )



ggsave(file.path("output", "vel_temp_combined.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)


# test similarity between cci and enso values


vt_combined_data_forcorr <- vel_temp %>%
  select(date, mean_vcur, mean_temp) %>%
  drop_na()

correlation_vt <- cor(vt_combined_data_forcorr$mean_vcur, vt_combined_data_forcorr$mean_temp, method = "pearson")
print(paste("Pearson correlation between Mean Temp and Mean Velocity:", correlation_vt))

# anova

anova_result_vt <- aov(mean_vcur ~ mean_temp, data = vel_temp)

summary(anova_result_vt)




# linear regression of EAC CCI against SOI

# Fit a linear model
model_vt <- lm(mean_vcur ~ mean_temp, data = vel_temp)

# Summary of the model
summary(model_vt)



# Save diagnostic plots to a PNG file
png(filename = "output/vel-temp-diagnostic.png", width = 1200, height = 600)

# Set up 2x2 layout and plot diagnostics
par(mfrow = c(2, 2))
plot(model_vt)

# Close the PNG device
dev.off()


# Cross-correlation function to explore lagged relationships
ccf_result <- ccf(vel_temp$mean_vcur, vel_temp$mean_temp, lag.max = 12, plot = TRUE,
                  main = "Cross-Correlation between Mean Monthly velocity and Mean Monthly Temperature")



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
vel_temp$temp_lagged <- dplyr::lag(vel_temp$mean_temp, n = abs(best_lag))

# If lag is negative, shift EAC_CCI instead
if (best_lag < 0) {
  vel_temp$vel_lagged <- dplyr::lag(vel_temp$mean_vcur, n = abs(best_lag))
  lag_model_vt <- lm(vel_lagged ~ mean_temp, data = vel_temp)
} else {
  lag_model_vt <- lm(mean_vcur ~ temp_lagged, data = vel_temp)
}

# View regression summary
summary(lag_model_vt)





# Compare strength with index ---------------------------------------------



# plot with the index values

with_vel <- month_data_t %>%
  mutate(date = as.Date(trip_month)) %>%
  left_join(data_tbl, by = "date")


velocity_scale <- max(with_vel$eac_cci, na.rm = TRUE) / max(abs(with_vel$mean_vcur), na.rm = TRUE)

with_vel <- with_vel %>%
  mutate(vel_scale = mean_vcur * velocity_scale)

p <- ggplot(with_vel, aes(x = date)) +
  geom_line(aes(y = eac_cci), colour = "black", linewidth = 1) +
  geom_point(aes(y = eac_cci), size = 3, shape = 21, fill = "black") +
  geom_smooth(aes(y = eac_cci), method = "loess", se = TRUE, color = "black", linetype = "dashed") +
  geom_line(aes(y = mean_vcur), colour = "blue", linewidth = 1) +
  geom_point(aes(y = mean_vcur), size = 3, shape = 21, fill = "blue") +
  geom_smooth(aes(y = mean_vcur), method = "loess", se = TRUE, color = "blue", linetype = "dashed") +
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



ggsave(file.path("output", "eac_cci_vel_combined.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)
