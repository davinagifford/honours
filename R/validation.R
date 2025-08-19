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

# Plot
ggplot(data_tbl, aes(x = date, y = mean_vcur)) +
  geom_line(color = "black") +
  geom_point(color = "black") +
  labs(
    title = "Monthly Mean VCUR",
    x = "Date",
    y = "Mean VCUR"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




# plot the data
p <- ggplot(data_tbl, aes(x = date, y = mean_vcur)) +
  geom_line(linewidth = 1, colour = "black") +
  geom_point(size = 3, shape = 21, fill = "black") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    x = "Time",
    y = expression("Monthly Mean Velocity (m/s"^{-1}*")"),
    title = "Mean Velocity Over Time"
  ) +
  theme(axis.title = element_text(size = 14), 
        plot.title = element_text(size = 18, face = "bold"))

ggsave(file.path("output", "mean-velocity_100m.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)


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
                     sec.axis = sec_axis(~ . / velocity_scale,  name = "Mean velocity")) +
  theme(
    axis.title.y = element_text(size = 14, color = "black"),
    axis.title.y.right = element_text(size = 14, color = "blue"),
    plot.title = element_text(size = 18, face = "bold")
  ) +
  labs(
    x = "Time",
    title = "EAC copepod composition index & Monthly mean velocity ")



ggsave(file.path("output", "eac_cci_vel_combined.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)
