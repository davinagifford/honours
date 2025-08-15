### eac_cci_clim_comparison.R
###
### plot the various climate indicators for comparison with the EAC CCI
###
### Created: 2025-08-14
### Author: Davina Gifford
### Last updated: 2025-08-14
### Edited by: Davina Gifford

library(tidyverse)
library(ggplot2)
library(lubridate)

# load SOI data

soi_data <- read_csv("data/soi_monthly.csv")

# Convert yyyymm to Date (first day of month)
soi_data <- soi_data %>%
  mutate(Date = ymd(paste0(Date, "01")))

p <- ggplot(soi_data, aes(x = Date, y = SOI)) +
  geom_line(linewidth = 1, colour = "black") +
  geom_point(size = 3, shape = 21, fill = "black") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  geom_smooth(method = "loess", se = TRUE, color = "black", linetype = "dashed") +
  theme(axis.title = element_text(size = 14), 
        plot.title = element_text(size = 18, face = "bold")) +
  labs(x = "Time",
       y = "SOI",
       title = "Southern Oscillation Index ")

ggsave(file.path("output", "SOI.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)


# load SAM data

sam_data <- read_csv("data/AAO-SAM.csv")

sam_data <- sam_data %>%
  filter(year > 2009)

# Calculate monthly average
sam_monthly <- sam_data %>%
  group_by(year, month) %>%
  summarise(
    mean_sam = mean(aao_index_cdas, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    full_date = as.Date(paste(year, month, "01", sep = "-"))
  )

p <- ggplot(sam_monthly, aes(x = full_date, y = mean_sam)) +
  geom_line(linewidth = 1, colour = "black") +
  geom_point(size = 3, shape = 21, fill = "black") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  geom_smooth(method = "loess", se = TRUE, color = "black", linetype = "dashed") +
  theme(axis.title = element_text(size = 14), 
        plot.title = element_text(size = 18, face = "bold")) +
  labs(x = "Time",
       y = "SAM (Monthly Mean)",
       title = "Southern Annular Mode (SAM) - Monthly Average")

ggsave(file.path("output", "SAM.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)


# ENSO

 enso_data <- read_csv("data/nina34.anom.csv")

 enso_data <- enso_data %>% 
   filter(Date > "2010-01-01") %>% # Filter for dates after 2009
   filter(Date < "2025-01-01") %>% # Filter for dates before 2025
   mutate(anom = `Nino Anom 3.4 Index  using ersstv5 from CPC  missing value -99.99 https://psl.noaa.gov/data/timeseries/month/`) 
 
 head(enso_data)
 p <- ggplot(enso_data, aes(x = Date, y = anom )) +
   geom_line(linewidth = 1, colour = "black") +
   geom_point(size = 3, shape = 21, fill = "black") +
   scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
   geom_smooth(method = "loess", se = TRUE, color = "black", linetype = "dashed") +
   theme(axis.title = element_text(size = 14), 
         plot.title = element_text(size = 18, face = "bold")) +
   labs(x = "Time",
        y = "ENSO Index",
        title = "ENSO Index (NINO3.4 Anomaly)") 
 ggsave(file.path("output", "ENSO.png"),
        plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
        device = png)
 