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

# calculate a climatology


climatology_str <- data_tbl %>%
  group_by(month) %>%
  summarise(
    str_clim = mean(mean_vcur, na.rm = TRUE), .groups = "drop"
  )



# Strength against climatology --------------------------------------------

strength_with_clim <- data_tbl %>%
  left_join(climatology_str, by = "month") %>%
  mutate(
    vcur_anomaly = mean_vcur - str_clim,
    anom_label = case_when(
      vcur_anomaly > 0 ~ "Anomaly (+)",
      TRUE ~ "Anomaly (-)"
    )
  )

str_anom <- strength_with_clim %>% 
  select(vcur_anomaly)

# Plot strength with climatology
p <- ggplot(strength_with_clim, aes(x = date, y = mean_vcur)) +
  geom_segment(aes(xend = date, yend = str_clim,
                   colour = anom_label, linetype = anom_label),
               linewidth = 1.5) +
  geom_line(aes(y = str_clim, colour = "Climatology", linetype = "Climatology"),
            linewidth = 1) +
  geom_point(size = 3, shape = 21, fill = "black") +
  scale_colour_manual(breaks = c("Anomaly (+)", "Anomaly (-)", "Climatology"),
                      values = c("red", "blue", "black")) +
  scale_linetype_manual(breaks = c("Anomaly (+)", "Anomaly (-)", "Climatology"),
                        values = c("solid", "solid", "solid")) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  theme(axis.title = element_text(size = 14),
        plot.title = element_text(size = 18, face = "bold")) +
  labs(x = "Time",
       y = "Strength ",
       colour = NULL,
       linetype = NULL,
       title = "EAC strength (Monthly Average)")

  
ggsave(file.path("output", "strength-with-clim.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)



# # is there a relationship between velocity and temperature?  ------------



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

with_vel <- anomaly %>%
  mutate(date = as.Date(trip_month)) %>%
  full_join(data_tbl, by = "date")


velocity_scale <- max(with_vel$anomaly, na.rm = TRUE) / max(abs(with_vel$mean_vcur), na.rm = TRUE)

with_vel <- with_vel %>%
  mutate(vel_scale = mean_vcur * velocity_scale)

p <- ggplot(with_vel, aes(x = date)) +
  geom_line(aes(y = anomaly), colour = "black", linewidth = 1) +
  geom_point(aes(y = anomaly), size = 3, shape = 21, fill = "black") +
  geom_smooth(aes(y = anomaly), method = "loess", se = TRUE, color = "black", linetype = "dashed") +
  geom_line(aes(y = mean_vcur), colour = "blue", linewidth = 1) +
  geom_point(aes(y = mean_vcur), size = 3, shape = 21, fill = "blue") +
  geom_smooth(aes(y = mean_vcur), method = "loess", se = TRUE, color = "blue", linetype = "dashed") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(name = "EAC copepod composition index - anomaly",
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

# test similarity between cci and strength values

with_vel <- with_vel %>% 
  mutate(
    year = year(trip_month),
    month = month(trip_month)
  ) 

str_combined_data_forcorr <- with_vel %>%
  select(month, observed_eac_cci, mean_vcur) %>%
  drop_na()

correlation_str <- cor(str_combined_data_forcorr$observed_eac_cci, str_combined_data_forcorr$mean_vcur, method = "pearson")
print(paste("Pearson correlation between EAC CCI and strength:", correlation_str)) 

# anova

anova_result_str <- aov(observed_eac_cci ~ mean_vcur, data = str_combined_data_forcorr)

summary(anova_result_str)

# linear regression of EAC CCI against strength

# Fit a linear model
model_str <- lm(observed_eac_cci ~ mean_vcur, data = with_vel)

# Summary of the model
summary(model_str)

anova(model_str)

# Plot the regression
plot(with_vel$mean_vcur, with_vel$observed_eac_cci, main = "Regression of EAC CCI on strength",
     xlab = "Mean Strength", ylab = "EAC CCI", pch = 19)
abline(model_str, col = "blue", lwd = 2)


plot(model_str$fitted.values, model_str$residuals,
     main = "Residuals vs Fitted Values",
     xlab = "Fitted Values",
     ylab = "Residuals",
     pch = 19, col = "darkgreen")
abline(h = 0, col = "red", lwd = 2)


# Save diagnostic plots to a PNG file
png(filename = "output/eac-str-diagnostic.png", width = 1200, height = 600)

# Set up 2x2 layout and plot diagnostics
par(mfrow = c(2, 2))
plot(model_str)

# Close the PNG device
dev.off()


# Remove rows with NA in either column
clean_data <- na.omit(with_vel[, c("mean_vcur", "observed_eac_cci")])

# Run CCF on cleaned data
ccf_result_str <- ccf(clean_data$mean_vcur, clean_data$observed_eac_cci, lag.max = 12, plot = TRUE,
                      main = "Cross-Correlation between Strength and EAC CCI")



# Extract correlation values and corresponding lags
correlations <- ccf_result_str$acf[,1,1]
lags <- ccf_result_str$lag[,1,1]

# Find the lag with the highest absolute correlation
max_index <- which.max(abs(correlations))
best_lag <- lags[max_index]
best_corr <- correlations[max_index]

# Print result
cat("Lag with highest cross-correlation:", best_lag, "months\n")
cat("Correlation coefficient:", round(best_corr, 3), "\n")




# Create lagged strength variable
with_vel$str_lagged <- dplyr::lag(with_vel$mean_vcur, n = abs(best_lag))

# If lag is negative, shift EAC_CCI instead
if (best_lag < 0) {
  with_vel$eac_cci_lagged <- dplyr::lag(with_vel$observed_eac_cci, n = abs(best_lag))
  lag_model_str <- lm(eac_cci_lagged ~ mean_vcur, data = with_vel)
} else {
  lag_model_str <- lm(observed_eac_cci ~ str_lagged, data = with_vel)
}

# View regression summary
summary(lag_model_str)


cor.test(with_vel$observed_eac_cci, with_vel$str_lagged, method = "pearson", use = "complete.obs")



# compare anomaly of index with anomaly of strength -----------------------

cci_strn <- anomaly %>%
  select(date = trip_month, cci_anom = anomaly) %>%
  left_join(str_anom %>% mutate(date = data_tbl$date), by = "date")

anom_combined_data_forcorr <- cci_strn %>%
  select(date, cci_anom, vcur_anomaly) %>%
  drop_na()

correlation_anom <- cor(anom_combined_data_forcorr$cci_anom, anom_combined_data_forcorr$vcur_anomaly, method = "pearson")
print(paste("Pearson correlation between EAC CCI and strength anomalies:", correlation_anom)) 

# anova

anova_result_anom <- aov(cci_anom ~ vcur_anomaly, data = anom_combined_data_forcorr)

summary(anova_result_anom)

# linear regression of EAC CCI against strength anomalies

# Fit a linear model
model_anom <- lm(cci_anom ~ vcur_anomaly, data = anom_combined_data_forcorr)

# Summary of the model
summary(model_anom)

# Plot the regression
plot(anom_combined_data_forcorr$vcur_anomaly, anom_combined_data_forcorr$cci_anom, main = "Regression of EAC CCI anomalies on strength anomalies",
     xlab = "Mean Strength anomalies", ylab = "EAC CCI anomalies", pch = 19)
abline(model_anom, col = "blue", lwd = 2)


plot(model_anom$fitted.values, model_anom$residuals,
     main = "Residuals vs Fitted Values",
     xlab = "Fitted Values",
     ylab = "Residuals",
     pch = 19, col = "darkgreen")
abline(h = 0, col = "red", lwd = 2)


# Save diagnostic plots to a PNG file
png(filename = "output/anom-diagnostic.png", width = 1200, height = 600)

# Set up 2x2 layout and plot diagnostics
par(mfrow = c(2, 2))
plot(model_anom)

# Close the PNG device
dev.off()


# Remove rows with NA in either column
clean_data <- na.omit(with_vel[, c("mean_vcur", "observed_eac_cci")])

# Run CCF on cleaned data
ccf_result_anom <- ccf(anom_combined_data_forcorr$vcur_anomaly, anom_combined_data_forcorr$cci_anom, lag.max = 12, plot = TRUE,
                      main = "Cross-Correlation between the anomalies in Strength and EAC CCI")



# Extract correlation values and corresponding lags
correlations <- ccf_result_anom$acf[,1,1]
lags <- ccf_result_anom$lag[,1,1]

# Find the lag with the highest absolute correlation
max_index <- which.max(abs(correlations))
best_lag <- lags[max_index]
best_corr <- correlations[max_index]

# Print result
cat("Lag with highest cross-correlation:", best_lag, "months\n")
cat("Correlation coefficient:", round(best_corr, 3), "\n")




# Create lagged strength variable
anom_combined_data_forcorr$vcur_lagged <- dplyr::lag(anom_combined_data_forcorr$vcur_anomaly, n = abs(best_lag))

# If lag is negative, shift EAC_CCI instead
if (best_lag < 0) {
  anom_combined_data_forcorr$cci_anom_lagged <- dplyr::lag(anom_combined_data_forcorr$cci_anom, n = abs(best_lag))
  lag_model_anom <- lm(cci_anom_lagged ~ vcur_anomaly, data = anom_combined_data_forcorr)
} else {
  lag_model_anom <- lm(cci_anom ~ vcur_lagged, data = anom_combined_data_forcorr)
}

# View regression summary
summary(lag_model_anom)


ggplot(anom_combined_data_forcorr, aes(x = vcur_anomaly, y = cci_anom)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(
    x = "EAC Strength Anomaly (VCUR)",
    y = "EAC CCI Anomaly",
    title = "Relationship Between EAC Strength and CCI Anomalies"
  )

# compare temp with sst ---------------------------------------------------

# read sst data
sst <- samples %>% 
  select(sample_time, sst)
colnames(sst) <- c("date", "sst")
sst <- sst %>%
    mutate(
      date = as.Date(sub("T.*", "", date)),
      year = year(date),
      month = month(date)
    ) %>%
      group_by(year, month) %>%
      summarise(sst = mean(sst, na.rm = TRUE), .groups = "drop") %>% 
      mutate(date = as.Date(paste(year, month, "01", sep = "-"))) 


sst_temp <- sst %>% 
  left_join(vel_temp %>% mutate(date = vel_temp$date), by = "date")

sst_temp <- sst_temp %>% 
  select(year, month, date, sst, mean_temp)


correlation_sst <- cor(sst_temp$sst, sst_temp$mean_temp, method = "pearson", use = "complete.obs")
print(paste("Pearson correlation between Index SST and observed Temperature:", correlation_sst)) 

# anova

anova_result_sst <- aov(sst ~ mean_temp, data = sst_temp)

summary(anova_result_sst)

# linear regression of EAC CCI against strength anomalies

# Fit a linear model
model_sst <- lm(sst ~ mean_temp, data = sst_temp)

# Summary of the model
summary(model_sst)

# Plot the regression
plot(anom_combined_data_forcorr$vcur_anomaly, anom_combined_data_forcorr$cci_anom, main = "Regression of EAC CCI anomalies on strength anomalies",
     xlab = "Mean Strength anomalies", ylab = "EAC CCI anomalies", pch = 19)
abline(model_anom, col = "blue", lwd = 2)


plot(model_anom$fitted.values, model_anom$residuals,
     main = "Residuals vs Fitted Values",
     xlab = "Fitted Values",
     ylab = "Residuals",
     pch = 19, col = "darkgreen")
abline(h = 0, col = "red", lwd = 2)


# Save diagnostic plots to a PNG file
png(filename = "output/anom-diagnostic.png", width = 1200, height = 600)

# Set up 2x2 layout and plot diagnostics
par(mfrow = c(2, 2))
plot(model_anom)

# Close the PNG device
dev.off()


# Remove rows with NA in either column
clean_data <- na.omit(with_vel[, c("mean_vcur", "observed_eac_cci")])

# Run CCF on cleaned data
ccf_result_anom <- ccf(anom_combined_data_forcorr$vcur_anomaly, anom_combined_data_forcorr$cci_anom, lag.max = 12, plot = TRUE,
                       main = "Cross-Correlation between the anomalies in Strength and EAC CCI")



# Extract correlation values and corresponding lags
correlations <- ccf_result_anom$acf[,1,1]
lags <- ccf_result_anom$lag[,1,1]

# Find the lag with the highest absolute correlation
max_index <- which.max(abs(correlations))
best_lag <- lags[max_index]
best_corr <- correlations[max_index]

# Print result
cat("Lag with highest cross-correlation:", best_lag, "months\n")
cat("Correlation coefficient:", round(best_corr, 3), "\n")




# Create lagged strength variable
anom_combined_data_forcorr$vcur_lagged <- dplyr::lag(anom_combined_data_forcorr$vcur_anomaly, n = abs(best_lag))

# If lag is negative, shift EAC_CCI instead
if (best_lag < 0) {
  anom_combined_data_forcorr$cci_anom_lagged <- dplyr::lag(anom_combined_data_forcorr$cci_anom, n = abs(best_lag))
  lag_model_anom <- lm(cci_anom_lagged ~ vcur_anomaly, data = anom_combined_data_forcorr)
} else {
  lag_model_anom <- lm(cci_anom ~ vcur_lagged, data = anom_combined_data_forcorr)
}

# View regression summary
summary(lag_model_anom)




# compare SST with Strength -----------------------------------------------

sst_str <- sst %>% 
  left_join(vel_temp %>% mutate(date = vel_temp$date), by = "date")

sst_str <- sst_str %>% 
  select(year, month, date, sst, mean_vcur)

correlation_sst_str <- cor(sst_str$sst, sst_str$mean_vcur, method = "pearson", use = "complete.obs")
print(paste("Pearson correlation between Index SST and observed strength:", correlation_sst_str)) 

# anova

anova_result_sst_str <- aov(sst ~ mean_vcur, data = sst_str)

summary(anova_result_sst_str)

# linear regression of EAC CCI against strength anomalies

# Fit a linear model
model_sst_str <- lm(sst ~ mean_vcur, data = sst_str)

# Summary of the model
summary(model_sst_str)


# Remove rows with NA in either column
clean_data_sst_str <- na.omit(sst_str[, c("mean_vcur", "sst")])

# Run CCF on cleaned data
ccf_result_sst_str <- ccf(clean_data_sst_str$mean_vcur, clean_data_sst_str$sst, lag.max = 12, plot = TRUE,
                       main = "Cross-Correlation between the Mean strength and SST used in Index")



# Extract correlation values and corresponding lags
correlations <- ccf_result_sst_str$acf[,1,1]
lags <- ccf_result_sst_str$lag[,1,1]

# Find the lag with the highest absolute correlation
max_index <- which.max(abs(correlations))
best_lag <- lags[max_index]
best_corr <- correlations[max_index]

# Print result
cat("Lag with highest cross-correlation:", best_lag, "months\n")
cat("Correlation coefficient:", round(best_corr, 3), "\n")




# Create lagged strength variable
sst_str$vcur_lagged <- dplyr::lag(sst_str$mean_vcur, n = abs(best_lag))

# If lag is negative, shift SST instead
if (best_lag < 0) {
  sst_str$sst_lagged <- dplyr::lag(sst_str$sst, n = abs(best_lag))
  lag_model_sst_str <- lm(sst_lagged ~ mean_vcur, data = sst_str)
} else {
  lag_model_sst_str <- lm(sst ~ vcur_lagged, data = sst_str)
}

# View regression summary
summary(lag_model_sst_str)

cor(sst_str$sst_lagged, sst_str$mean_vcur, method = "pearson", use = "complete.obs")
cor.test(sst_str$sst_lagged, sst_str$mean_vcur, method = "pearson", use = "complete.obs")


# compare all the things --------------------------------------------------

# Combine temp, sst, strength, and index data

all_data <- with_vel %>%
  select(date, observed_eac_cci, mean_vcur) %>%
  left_join(vel_temp %>% select(date, mean_temp), by = "date") %>%
  left_join(sst %>% select(date, sst), by = "date")

# create a scatterplot matrix of the variables

pairs(all_data[, -1], main = "Scatterplot Matrix of EAC CCI, Strength, Temperature, and SST")
# calculate correlation matrix
cor_matrix <- cor(all_data[, -1], use = "pairwise.complete.obs")
print(cor_matrix)
# visualize correlation matrix
library(corrplot)
corrplot(cor_matrix, method = "circle", type = "upper", tl.col = "black", tl.srt = 45,
         title = "Correlation Matrix of EAC CCI, Strength, Temperature, and SST", mar = c(0,0,1,0))


cor(all_data$sst, all_data$observed_eac_cci, method = "pearson", use = "complete.obs")

summary(lm(sst ~ observed_eac_cci, data = all_data))
summary(lm(sst ~ mean_temp, data = all_data))
summary(lm(sst ~ mean_vcur, data = all_data))
summary(lm(mean_temp ~ observed_eac_cci, data = all_data))
summary(lm(mean_temp ~ mean_vcur, data = all_data))
summary(lm(mean_vcur ~ observed_eac_cci, data = all_data))
summary(lm(observed_eac_cci ~ mean_vcur, data = all_data))
summary(lm(observed_eac_cci ~ mean_temp, data = all_data))
summary(lm(observed_eac_cci ~ sst, data = all_data))
summary(lm(mean_vcur ~ sst, data = all_data))
summary(lm(mean_temp ~ sst, data = all_data))
summary(lm(mean_vcur ~ mean_temp, data = all_data))
summary(lm(mean_vcur ~ mean_temp, data = all_data))
summary(lm(north_index ~ mean_vcur, data = all_data))






# Validate the strength against just the index values from the North

all_data <- all_data %>% 
  left_join(month_data_north %>% select(date = trip_month, north_index = eac_cci), by = "date")


pairs(all_data[, -1], main = "Scatterplot Matrix of EAC CCI, Strength, Temperature, and SST")
# calculate correlation matrix
cor_matrix <- cor(all_data[, -1], use = "pairwise.complete.obs")
print(cor_matrix)


# Run Pearson correlation
correlation <- cor(all_data$mean_vcur, all_data$north_index, method = "pearson", use = "pairwise.complete.obs")
cor_test <- cor.test(all_data$mean_vcur, all_data$north_index, method = "pearson", use = "pairwise.complete.obs")

# Output results
print(paste("Pearson correlation:", round(correlation, 3)))
print(cor_test)
