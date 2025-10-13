### compare_27.R

### Compare the Index, observation data and Copernicus data from 27S
###
### Created: 2025-09-15
### Author: Davina Gifford
### Last updated: 2025-09-15
### Edited by: Davina Gifford


# load libraries

library(tidyverse)
library(ggplot2)
library(corrplot)
library(mgcv)


# filter copernicus model data to 27 s
# use model data from validation_cop_model.R

cop_27 <- bind_rows(model_vel_full_50, raw_model_vel)

cop_27$depth <- as.numeric(cop_27$depth)


cop_27 <- cop_27 %>%
  filter(latitude == -27) %>% 
  filter(depth < 60) %>%
  filter(longitude > 153.9) %>% 
  filter(longitude < 155.0) %>% 
  mutate(
    date = as.Date(sub("T.*", "", time)),
    year = year(date),
    month = month(date)
  ) %>%
  group_by(year, month) %>%
  summarise(mean_vel = mean(vo, na.rm = TRUE), .groups = "drop") %>% 
  mutate(date = as.Date(paste(year, month, "01", sep = "-"))) 

# convert velocity to positive strength value
cop_27 <- cop_27 %>%
  mutate(mean_vel = abs(mean_vel))

# climatology at 27 for Copernicus data -----------------------------------

cop_clim_27 <-  cop_27 %>%
  mutate(
    doy = yday(date)
  )

lm_fit_cop_27 <- gam(mean_vel ~ s(doy, bs = "cc", k = 5),
                     data = cop_clim_27,
                     method = "REML",
                     knots = list(doy = c(0, 365)))

mean_time_cop_27 <- median(cop_clim_27$date)
time0_cop_27 <- min(cop_clim_27$date)

climatology_cop_27 <-
  tibble(date = floor_date(mean_time_cop_27, unit = "year") + days(0:364),
         doy = yday(date),
         time_x = time_length(interval(time0_cop_27, mean_time_cop_27), unit = "day"))

clim_terms_cop_27 <- predict(lm_fit_cop_27, type = "terms", newdata = climatology_cop_27)

head(clim_terms_cop_27)

climatology_cop_27 <-
  climatology_cop_27 %>%
  mutate(doy_eff = clim_terms_cop_27[, "s(doy)"],
         intercept = coef(lm_fit_cop_27)["(Intercept)"],
         str_clim = doy_eff + intercept) %>%
  select(!c(doy_eff, intercept))

climatology_cop_27 <- climatology_cop_27 %>%
  mutate(date = as.Date(date))


month_clim_cop_27 <- climatology_cop_27 %>%
  mutate(month = month(date)) %>%
  group_by(month) %>%
  summarise(month_clim_str = mean(str_clim, na.rm = TRUE), .groups = "drop") %>% 
  rename(month_clim_cop_27 = month_clim_str) 


# add in CCI climatology
clim_compare_short_cop_27 <- month_clim_cop_27 %>%
  select(month, month_clim_cop_27) %>%
  left_join(month_clim_str, by = "month") %>% # Observed data
  left_join(month_climatology, by = "month") # CCI






# Pivot longer for easier plotting
clim_compare_long_cop_27 <- clim_compare_short_cop_27 %>%
  pivot_longer(cols = c(month_clim_cop_27, month_clim_str, month_clim),
               names_to = "climatology_type",
               values_to = "value")

# Calculate correlations

# correlation between observation data at 27 and copernicus data at 27

# Calculate correlation


# Run Pearson correlation
correlation_27 <- cor(clim_compare_short_cop_27$month_clim_cop_27, clim_compare_short_cop_27$month_clim_str, method = "pearson", use = "complete.obs")
cor_test_27 <- cor.test(clim_compare_short_cop_27$month_clim_cop_27, clim_compare_short_cop_27$month_clim_str, method = "pearson", use = "complete.obs")

# Output results
print(paste("Pearson correlation:", round(correlation_27, 3)))
print(cor_test_27)


p_val_27 <- cor_test_27$p.value



cor_matrix_cop_27 <- clim_compare_short_cop_27 %>%
  select(month_clim_cop_27, month_clim_str, month_clim) %>%
  cor(method = "pearson")



print(cor_matrix_cop_27)


# Plot the correlation matrix
corrplot(cor_matrix_cop_27, method = "circle", type = "upper", tl.cex = 0.8)


# Plot comparison
p <- ggplot(clim_compare_long_cop_27, aes(x = month, y = value, color = climatology_type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3, shape = 21, fill = "white") +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_color_manual(values = c("month_clim_str" = "blue", "month_clim" = "red", "month_clim_cop_27" = "darkgreen"),
                     labels = c("EAC CCI Climatology", "Modelled Strength Climatology", "Observed Strength Climatology")) +
  labs(
    x = "Month",
    y = "Climatology Value",
    color = "Climatology Type",
    title = "Comparison of EAC Strength (observed and modelled) and EAC CCI Climatologies at 27S",
    #subtitle = paste("Pearson correlation:", round(correlation_sub_k, 3), "| p-value:", signif(p_val_sub_k, 3))
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )


ggsave(file.path("output", "climatology-comparison_cop_27.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)


# plot a scatterplot of the two

ggplot(data = clim_compare_short_sub) +
  geom_point(mapping = aes(x = month_clim, y = month_clim_str)) +
  geom_smooth(mapping = aes(x = month_clim, y = month_clim_str), method = "lm", se = TRUE, color = "blue", linetype = "dashed") +
  labs(
    x = "EAC CCI Climatology",
    y = "Current Strength Climatology",
    title = "Scatterplot of Current Strength vs EAC CCI Climatologies",
    subtitle = paste("Pearson correlation:", round(correlation_sub, 3), "| p-value:", signif(p_val_sub, 3))
  ) 

