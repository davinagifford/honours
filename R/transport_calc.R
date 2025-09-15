library(ncdf4)
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyverse)
fn <- "D:/HONOURS/DavinaG_2025_Honours/data/cmems_mod_glo_phy_my_0.083deg_P1M-m_1757126612458.nc"

# ---------- helpers ----------
cell_thickness_from_centres <- function(z) {
  n <- length(z); stopifnot(n >= 2)
  dzm <- diff(z)
  dz  <- c(dzm[1], (dzm[-1] + dzm[-length(dzm)])/2, dzm[length(dzm)])
  dz
}
deg_lon_to_m <- function(dlon, lat_deg) 111320 * cos(lat_deg*pi/180) * dlon
nearest_idx <- function(v, target) which.min(abs(v - target))

# ---------- coords & time ----------
nc <- nc_open(fn)
lon <- nc$dim$longitude$vals
lat <- nc$dim$latitude$vals
dep <- nc$dim$depth$vals
t_raw <- nc$dim$time$vals
t_units <- nc$dim$time$units

origin_str <- sub(".*since\\s+", "", t_units)
origin <- suppressWarnings(ymd_hms(origin_str, quiet=TRUE)); if (is.na(origin)) origin <- suppressWarnings(ymd(origin_str, quiet=TRUE))
sec_per <- if (grepl("^days", t_units)) 86400 else if (grepl("^hours", t_units)) 3600 else 1
time <- origin + t_raw * sec_per

# Dimension order of 'vo' in THIS file
vo_info <- nc$var$vo
dim_names <- vapply(vo_info$dim, function(d) d$name, character(1))
dim_lengths <- vapply(vo_info$dim, function(d) d$len, integer(1))
# Typically: c("time","depth","latitude","longitude") for CMEMS monthly, but we won't assume

# ---------- selection ----------
lon_min <- 153; lon_max <- 155
z_max <- 60 # values less than 60 m depth  
southward_positive <- TRUE

idx_lon <- which(lon >= lon_min & lon <= lon_max)
idx_dep <- which(dep >= 0 & dep <= z_max)

# targets
targets <- c(`27S`=-27, `32S`=-32)
res_list <- list()

for (nm in names(targets)) {
  target_lat <- targets[[nm]]
  jlat <- nearest_idx(lat, target_lat)
  lat_used <- lat[jlat]
  
  # Build start/count vectors in the EXACT order of the 'vo' dimensions
  starts <- counts <- integer(length(dim_names))
  for (i in seq_along(dim_names)) {
    dn <- dim_names[i]
    if (dn == "time")   { starts[i] <- 1; counts[i] <- length(time) }
    if (dn == "depth")  { starts[i] <- min(idx_dep); counts[i] <- length(idx_dep) }
    if (dn == "latitude"){ starts[i] <- jlat; counts[i] <- 1 }
    if (dn == "longitude"){ starts[i] <- min(idx_lon); counts[i] <- length(idx_lon) }
  }
  
  # Pull subset (collapse_degen=FALSE keeps the singleton latitude dim)
  vo_sub <- ncvar_get(nc, "vo", start = starts, count = counts, collapse_degen = TRUE)
  
  # Figure which axis is which index in the retrieved array
  arr_dim <- dim(vo_sub)
  i_time <- which(arr_dim == length(time))[1]
  i_lon  <- which(arr_dim == length(idx_lon))[1]
  i_dep  <- which(arr_dim == length(idx_dep))[1]
  
  
  # Reorder to [time, lon, depth]
  vo_tld <- aperm(vo_sub, perm = c(i_time, i_lon, i_dep))
  
  # Prepare metrics
  lons_use <- lon[idx_lon]
  deps_use <- dep[idx_dep]
  # dx from lon degrees at this latitude (same for all depths)
  dlon <- diff(lons_use)
  dx <- c(dlon[1], (dlon[-1] + dlon[-length(dlon)])/2, dlon[length(dlon)])
  dx_m <- deg_lon_to_m(dx, lat_used)
  dz_m <- cell_thickness_from_centres(deps_use)
  A <- outer(dx_m, dz_m)  # [lon,depth] m^2
  
  # Integrate each month
  nT <- dim(vo_tld)[1]
  Ts <- numeric(nT)
  for (it in seq_len(nT)) {
    v <- vo_tld[it,,]           # [lon,depth], northward m/s
    if (southward_positive) v <- -v
    Ts[it] <- sum(v * A, na.rm = TRUE) / 1e6  # Sv
  }
  
  res_list[[nm]] <- tibble(time = time, section = nm, latitude = lat_used, transport_Sv = Ts)
}

nc_close(nc)

transports <- bind_rows(res_list) |> arrange(time, section)

# Difference 27S - 32S
diff_df <- transports |>
  select(section, time, transport_Sv) |>
  tidyr::pivot_wider(names_from = section, values_from = transport_Sv) |>
  mutate(diff_27S_minus_32S = `27S` - `32S`)

# Quick plots
ggplot(transports, aes(time, transport_Sv, colour = section)) +
  geom_line() +
  labs(x = "Month", y = "Transport (Sv, poleward +)", colour = "Latitude",
       title = "EAC meridional transport (0–1500 m) from CMEMS reanalysis")

ggplot(diff_df, aes(time, diff_27S_minus_32S)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line() +
  labs(x = "Month", y = "Δ Transport (Sv)",
       title = "Transport change between 27°S and 32°S (0–1500 m)")



# check correlation between the two transports

cor_test_tran <- cor.test(
  transports$transport_Sv[transports$section == "27S"],
  transports$transport_Sv[transports$section == "32S"]
)
print(cor_test_tran)




# understand difference between two regions -------------------------------

library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)

# Wide form for paired comparisons
tw <- transports |>
  select(time, section, transport_Sv) |>
  pivot_wider(names_from = section, values_from = transport_Sv) |>
  arrange(time)

# Optional: drop months with missing at either latitude
tw <- tw |> filter(is.finite(`27S`) & is.finite(`32S`))

# Quick look
summary(tw)

# absolute downstream change in Sv

tw <- tw |>
  mutate(delta_Sv = `27S` - `32S`)


# Retention / percentage change

tw <- tw |>
  mutate(retention_pct = 100 * (`32S` / `27S`))


# seasonal structure vs anomalies

clim <- tw |>
  mutate(mon = month(time)) |>
  group_by(mon) |>
  summarise(T27_clim = mean(`27S`, na.rm=TRUE),
            T32_clim = mean(`32S`, na.rm=TRUE),
            dT_clim  = mean(delta_Sv, na.rm=TRUE),
            retention_clim = 100 * mean(`32S`/`27S`, na.rm=TRUE),
            .groups="drop")

# Plot climatology
ggplot(clim, aes(mon, dT_clim)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line() + geom_point() +
  labs(x = "Month", y = "Δ Transport climatology (Sv)",
       title = "Seasonal mean change: 27°S – 32°S")


#Anomalies (to study variability independent of the mean seasonal cycle)

tw_anom <- tw |>
  mutate(mon = month(time)) |>
  left_join(clim |> select(mon, T27_clim, T32_clim), by="mon") |>
  mutate(T27_anom = `27S` - T27_clim,
         T32_anom = `32S` - T32_clim,
         dT_anom  = T27_anom - T32_anom,
         ret_anom = 100 * ( (`32S`/`27S`) - (T32_clim/T27_clim) ))

# Correlation of anomalies between latitudes
cor_Tanom <- with(tw_anom, cor(T27_anom, T32_anom, use="pairwise"))
cor_Tanom


# Downstream relationship & lag
# (a) Instantaneous linear relation (how much of 27°S appears at 32°S)

fit_inst <- lm(`32S` ~ `27S`, data = tw)
summary(fit_inst)
# Slope ~ "retention factor" on average (not the same as mean ratio, but intuitive)

# (b) Propagation lag (does 32°S follow 27°S by N months?)
cc <- ccf(tw_anom$T27_anom, tw_anom$T32_anom, lag.max = 12, plot = FALSE)
lag_months <- cc$lag[which.max(cc$acf)]
lag_months  # e.g., 1–3 months is common for downstream advection signals

k <- as.integer(lag_months)
if (k > 0) {
  tw_lag <- tw_anom |>
    mutate(T27_lead = dplyr::lead(T27_anom, k)) |>
    filter(is.finite(T27_lead) & is.finite(T32_anom))
  fit_lag <- lm(T32_anom ~ T27_lead, data = tw_lag)
  summary(fit_lag)
}


# Trend: is the downstream change itself changing over time?

# Linear trend
fit_trend <- lm(delta_Sv ~ time, data = tw)
summary(fit_trend)

# Robust or flexible alternative (GAM)
library(mgcv)
fit_gam <- gam(delta_Sv ~ s(as.numeric(time), k = 10), data = tw, method = "REML")
plot(fit_gam, residuals = TRUE, pch = 16, cex = 0.6, main = "ΔT smooth (Sv)")


# visuals

# Time series of T at both lats
tw |>
  pivot_longer(c(`27S`,`32S`), names_to = "section", values_to = "Sv") |>
  ggplot(aes(time, Sv, colour = section)) +
  geom_line() +
  labs(y = "Transport (Sv, poleward +)", colour = "Latitude",
       title = "EAC transport at 27°S vs 32°S") +
  theme_minimal()

# ΔT over time
ggplot(tw, aes(time, delta_Sv)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line() +
  labs(y = "Δ Transport (Sv) = 27°S − 32°S",
       title = "Downstream change in transport") +
  theme_minimal()

# Retention %
ggplot(tw, aes(time, retention_pct)) +
  geom_hline(yintercept = 100, linetype = 2) +
  geom_line() +
  labs(y = "Retention (%) = 100 × (T32 / T27)",
       title = "Fraction of 27°S transport seen again at 32°S") +
  theme_minimal()


# climatology

library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)

# assume `transports` has columns: time, section ("27S"/"32S"), transport_Sv

clim <- transports |>
  mutate(month = month(time, label = TRUE)) |>
  group_by(section, month) |>
  summarise(clim_val = mean(transport_Sv, na.rm=TRUE), .groups="drop")

ggplot(clim, aes(x = month, y = clim_val, colour = section, group = section)) +
  geom_line() + geom_point() +
  labs(title = "Climatology of EAC Transport",
       y = "Transport (Sv, poleward +)", x = "Month") +
  theme_minimal()


# anomaly style

# compute monthly climatology
clim_month <- transports |>
  mutate(month = month(time)) |>
  group_by(section, month) |>
  summarise(clim_val = mean(transport_Sv, na.rm=TRUE), .groups="drop")

# join back & calculate anomaly
anom <- transports |>
  mutate(month = month(time)) |>
  left_join(clim_month, by=c("section","month")) |>
  mutate(anomaly = transport_Sv - clim_val)

# plot
ggplot(anom, aes(time, transport_Sv)) +
  geom_line(aes(colour = section)) +
  geom_line(aes(y = clim_val, group=section), colour="black") +
  geom_segment(aes(x=time, xend=time,
                   y=clim_val, yend=transport_Sv,
                   colour = ifelse(anomaly>0, "pos","neg"))) +
  scale_colour_manual(values=c("pos"="red","neg"="blue","27S"="red","32S"="blue")) +
  labs(y = "Transport (Sv)", x = "Time",
       title = "Monthly Transport with Climatology & Anomalies") +
  theme_minimal()



library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)

# ----- 1) Monthly climatology for each section -----
clim <- transports |>
  mutate(mon = month(time, label = TRUE, abbr = TRUE)) |>
  group_by(section, mon) |>
  summarise(clim_Sv = mean(transport_Sv, na.rm = TRUE), .groups = "drop")

# ----- 2) ΔT (27S − 32S) climatology -----
clim_wide <- clim |>
  pivot_wider(names_from = section, values_from = clim_Sv)

clim_diff <- clim_wide |>
  mutate(delta_Sv = `27S` - `32S`) |>
  select(mon, delta_Sv)

# ----- 3) Plots -----
p_top <- ggplot(clim, aes(mon, clim_Sv, colour = section, group = section)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_colour_manual(values = c(`27S` = "#d73027", `32S` = "#4575b4")) +
  labs(
    title = "EAC Transport Climatologies",
    subtitle = "Monthly means; poleward (southward) transport is positive",
    x = "Month", y = "Transport (Sv)", colour = "Latitude"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

p_bottom <- ggplot(clim_diff, aes(mon, delta_Sv, group = 1)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    title = "ΔT Climatology (27°S − 32°S)",
    subtitle = "Positive = decrease downstream; Negative = increase downstream",
    x = "Month", y = "Δ Transport (Sv)"
  ) +
  theme_minimal(base_size = 12)

# ----- 4) Combine panels -----
# If you have patchwork installed:
# install.packages("patchwork")  # once
library(patchwork)
combined <- p_top / p_bottom + plot_layout(heights = c(2, 1))
combined

# ----- (Optional) Save to file -----
ggsave("transport_climatology_and_delta.png", combined, width = 10, height = 7, dpi = 200)


library(dplyr)
library(lubridate)

yearly_means <- transports |>
  mutate(year = year(time)) |>
  group_by(section, year) |>
  summarise(mean_transport_Sv = mean(transport_Sv, na.rm = TRUE), .groups = "drop")

print(yearly_means)

yearly_wide <- yearly_means |>
  tidyr::pivot_wider(names_from = section, values_from = mean_transport_Sv)

yearly_wide <- yearly_wide |>
  mutate(delta_Sv = `27S` - `32S`,
         retention_pct = 100 * (`32S` / `27S`))

print(yearly_wide)


library(ggplot2)

# Yearly means at each latitude
ggplot(yearly_means, aes(year, mean_transport_Sv, colour = section)) +
  geom_line() + geom_point() +
  labs(y = "Annual Mean Transport (Sv, poleward +)",
       x = "Year", colour = "Latitude",
       title = "Yearly Mean EAC Transport") +
  theme_minimal()

# Yearly ΔT
ggplot(yearly_wide, aes(year, delta_Sv)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line() + geom_point() +
  labs(y = "Δ Transport (Sv) = 27°S − 32°S",
       x = "Year",
       title = "Yearly Difference in Transport") +
  theme_minimal()


library(dplyr)
library(lubridate)
library(tidyr)

# 1. Yearly mean transport at each section
yearly_means <- transports |>
  mutate(year = year(time)) |>
  group_by(section, year) |>
  summarise(mean_transport_Sv = mean(transport_Sv, na.rm = TRUE), .groups = "drop")

# 2. Long-term mean per section (baseline for anomalies)
lt_mean <- yearly_means |>
  group_by(section) |>
  summarise(longterm_mean = mean(mean_transport_Sv, na.rm = TRUE), .groups = "drop")

# 3. Yearly anomalies (difference from long-term mean)
yearly_anoms <- yearly_means |>
  left_join(lt_mean, by = "section") |>
  mutate(anomaly_Sv = mean_transport_Sv - longterm_mean)

# 4. Wide format for easy comparison (adds delta and retention too)
yearly_wide <- yearly_anoms |>
  select(section, year, mean_transport_Sv, anomaly_Sv) |>
  pivot_wider(names_from = section,
              values_from = c(mean_transport_Sv, anomaly_Sv))

yearly_wide <- yearly_wide |>
  mutate(delta_Sv        = mean_transport_Sv_27S - mean_transport_Sv_32S,
         delta_anomaly   = anomaly_Sv_27S - anomaly_Sv_32S,
         retention_pct   = 100 * (mean_transport_Sv_32S / mean_transport_Sv_27S))

# Results
print(yearly_anoms)  # tidy table: section × year × mean × anomaly
print(yearly_wide)   # wide table: both sites + delta + retention

library(ggplot2)

# Annual anomalies by site
ggplot(yearly_anoms, aes(year, anomaly_Sv, colour = section)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line() + geom_point() +
  labs(y = "Annual Transport Anomaly (Sv)", x = "Year", colour = "Latitude",
       title = "EAC Transport Anomalies (annual means)") +
  theme_minimal()

# Δ anomaly between sites
ggplot(yearly_wide, aes(year, delta_anomaly)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line() + geom_point() +
  labs(y = "Δ Annual Anomaly (Sv) = 27°S − 32°S",
       x = "Year",
       title = "Downstream Change in Anomalies") +
  theme_minimal()


# assuming you already have yearly_wide from the previous step
# with columns anomaly_Sv_27S and anomaly_Sv_32S

# Correlation test
cor_test <- cor.test(yearly_wide$anomaly_Sv_27S,
                     yearly_wide$anomaly_Sv_32S,
                     use = "pairwise.complete.obs",
                     method = "pearson")

print(cor_test)
# cor_test$estimate gives correlation coefficient (r)
# cor_test$p.value gives significance

# Scatterplot with regression line
library(ggplot2)

ggplot(yearly_wide, aes(x = anomaly_Sv_27S, y = anomaly_Sv_32S)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_point(size = 3, colour = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, colour = "black") +
  labs(
    x = "Annual Transport Anomaly at 27°S (Sv)",
    y = "Annual Transport Anomaly at 32°S (Sv)",
    title = "Correlation of Annual Transport Anomalies (27°S vs 32°S)",
    subtitle = paste0("r = ", round(cor_test$estimate, 2),
                      ", p = ", signif(cor_test$p.value, 3))
  ) +
  theme_minimal()



# climatology for transport -----------------------------------------------
# 27 degree
transport_data_27 <- transports %>%
  filter(section == "27S") %>%
  mutate(
    doy = yday(time)
  )

lm_fit_trans_27 <- gam(transport_Sv ~ s(doy, bs = "cc", k = 5),
                       data = transport_data_27,
                       method = "REML",
                       knots = list(doy = c(0, 365)))

mean_time <- median(transport_data_27$time)
time0 <- min(transport_data_27$time)

climatology_trans_27 <-
  tibble(date = floor_date(mean_time, unit = "year") + days(0:364),
         doy = yday(date),
         time_x = time_length(interval(time0, mean_time), unit = "day"))

clim_terms_trans_27 <- predict(lm_fit_trans_27, type = "terms", newdata = climatology_trans_27)

head(clim_terms_trans_27)

climatology_trans_27 <-
  climatology_trans_27 %>%
  mutate(doy_eff = clim_terms_trans_27[, "s(doy)"],
         intercept = coef(lm_fit_trans_27)["(Intercept)"],
         trans_clim = doy_eff + intercept) %>%
  select(!c(doy_eff, intercept))

climatology_trans_27 <- climatology_trans_27 %>%
  mutate(date = as.Date(date))


month_clim_trans_27 <- climatology_trans_27 %>%
  mutate(month = month(date)) %>%
  group_by(month) %>%
  summarise(month_clim_trans_27 = mean(trans_clim, na.rm = TRUE), .groups = "drop")

clim_compare_short_trans_27 <- month_clim_trans_27 %>%
  select(month, month_clim_trans_27) %>%
  left_join(
    month_climatology,
    by = "month"
  )



# Pivot longer for easier plotting
clim_compare_long_trans_27 <- clim_compare_short_trans_27 %>%
  pivot_longer(cols = c(month_clim_trans_27, month_clim),
               names_to = "climatology_type",
               values_to = "value")

# Calculate correlation


# Run Pearson correlation
correlation_trans_27 <- cor(clim_compare_short_trans_27$month_clim, clim_compare_short_trans_27$month_clim_trans_27, method = "pearson", use = "complete.obs")
cor_test_trans_27 <- cor.test(clim_compare_short_trans_27$month_clim, clim_compare_short_trans_27$month_clim_trans_27, method = "pearson", use = "complete.obs")

# Output results
print(paste("Pearson correlation:", round(correlation_trans_27, 3)))
print(cor_test_trans_27)


p_val_trans_27 <- cor_test_trans_27$p.value



# Plot comparison
p <- ggplot(clim_compare_long_trans_27, aes(x = month, y = value, color = climatology_type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3, shape = 21, fill = "white") +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_color_manual(values = c("month_clim_trans_27" = "blue", "month_clim" = "red"),
                     labels = c("EAC CCI Climatology", "Transport Climatology")) +
  labs(
    x = "Month",
    y = "Climatology Value",
    color = "Climatology Type",
    title = "Comparison of EAC Transport and EAC CCI Climatologies",
    subtitle = paste("Pearson correlation:", round(correlation_trans_27, 3), "| p-value:", signif(p_val_trans_27, 3))
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )


ggsave(file.path("output", "climatology-comparison_trans_27.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)



# plot a scatterplot of the two

ggplot(data = clim_compare_short_trans_27) +
  geom_point(mapping = aes(x = month_clim, y = month_clim_trans_27)) +
  geom_smooth(mapping = aes(x = month_clim, y = month_clim_trans_27), method = "lm", se = TRUE, color = "blue", linetype = "dashed") +
  labs(
    x = "EAC CCI Climatology",
    y = "Transport Climatology",
    title = "Scatterplot of Transport vs EAC CCI Climatologies",
    subtitle = paste("Pearson correlation:", round(correlation_trans_27, 3), "| p-value:", signif(p_val_trans_27, 3))
  ) 




# 32 degree

transport_data_32 <- transports %>%
  filter(section == "32S") %>%
  mutate(
    doy = yday(time)
  )

lm_fit_trans_32 <- gam(transport_Sv ~ s(doy, bs = "cc", k = 5),
                       data = transport_data_32,
                       method = "REML",
                       knots = list(doy = c(0, 365)))

mean_time <- median(transport_data_32$time)
time0 <- min(transport_data_32$time)

climatology_trans_32 <-
  tibble(date = floor_date(mean_time, unit = "year") + days(0:364),
         doy = yday(date),
         time_x = time_length(interval(time0, mean_time), unit = "day"))

clim_terms_trans_32 <- predict(lm_fit_trans_32, type = "terms", newdata = climatology_trans_32)

head(clim_terms_trans_32)

climatology_trans_32 <-
  climatology_trans_32 %>%
  mutate(doy_eff = clim_terms_trans_32[, "s(doy)"],
         intercept = coef(lm_fit_trans_32)["(Intercept)"],
         trans_clim = doy_eff + intercept) %>%
  select(!c(doy_eff, intercept))

climatology_trans_32 <- climatology_trans_32 %>%
  mutate(date = as.Date(date))


month_clim_trans_32 <- climatology_trans_32 %>%
  mutate(month = month(date)) %>%
  group_by(month) %>%
  summarise(month_clim_trans_32 = mean(trans_clim, na.rm = TRUE), .groups = "drop")

clim_compare_short_trans_32 <- month_clim_trans_32 %>%
  select(month, month_clim_trans_32) %>%
  left_join(
    month_climatology,
    by = "month"
  )



# Pivot longer for easier plotting
clim_compare_long_trans_32 <- clim_compare_short_trans_32 %>%
  pivot_longer(cols = c(month_clim_trans_32, month_clim),
               names_to = "climatology_type",
               values_to = "value")

# Calculate correlation


# Run Pearson correlation
correlation_trans_32 <- cor(clim_compare_short_trans_32$month_clim, clim_compare_short_trans_32$month_clim_trans_32, method = "pearson", use = "complete.obs")
cor_test_trans_32 <- cor.test(clim_compare_short_trans_32$month_clim, clim_compare_short_trans_32$month_clim_trans_32, method = "pearson", use = "complete.obs")

# Output results
print(paste("Pearson correlation:", round(correlation_trans_32, 3)))
print(cor_test_trans_32)


p_val_trans_32 <- cor_test_trans_32$p.value



# Plot comparison
p <- ggplot(clim_compare_long_trans_32, aes(x = month, y = value, color = climatology_type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3, shape = 21, fill = "white") +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_color_manual(values = c("month_clim_trans_32" = "blue", "month_clim" = "red"),
                     labels = c("EAC CCI Climatology", "Transport Climatology")) +
  labs(
    x = "Month",
    y = "Climatology Value",
    color = "Climatology Type",
    title = "Comparison of EAC Transport and EAC CCI Climatologies",
    subtitle = paste("Pearson correlation:", round(correlation_trans_32, 3), "| p-value:", signif(p_val_trans_32, 3))
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )


ggsave(file.path("output", "climatology-comparison_trans_32.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)





# plot a scatterplot of the two

ggplot(data = clim_compare_short_trans_32) +
  geom_point(mapping = aes(x = month_clim, y = month_clim_trans_32)) +
  geom_smooth(mapping = aes(x = month_clim, y = month_clim_trans_32), method = "lm", se = TRUE, color = "blue", linetype = "dashed") +
  labs(
    x = "EAC CCI Climatology",
    y = "Transport Climatology",
    title = "Scatterplot of Transport vs EAC CCI Climatologies",
    subtitle = paste("Pearson correlation:", round(correlation_trans_32, 3), "| p-value:", signif(p_val_trans_32, 3))
  ) 


# Residuals  --------------------------------------------------------------
# linear regression of the two climatologies with the transport data


# Fit a linear model
trans_27_model <- lm(month_clim ~ month_clim_trans_27, data = clim_compare_short_trans_27)

# Summary of the model
summary(trans_27_model)



residuals_trans_27 <- resid(trans_27_model)
predicted_trans_27 <- fitted(trans_27_model)



plot(clim_compare_short_trans_27$month_clim_trans_27, trans_27_model$residuals,
     main = "Residuals vs Fitted Values",
     xlab = "Fitted Values",
     ylab = "Residuals",
     pch = 19, col = "darkgreen")
abline(h = 0, col = "red", lwd = 2)

# 32 degrees

# Fit a linear model
trans_32_model <- lm(month_clim ~ month_clim_trans_32, data = clim_compare_short_trans_32)

# Summary of the model
summary(trans_32_model)



residuals_trans_32 <- resid(trans_32_model)
predicted_trans_32 <- fitted(trans_32_model)



plot(clim_compare_short_trans_32$month_clim_trans_32, trans_32_model$residuals,
     main = "Residuals vs Fitted Values",
     xlab = "Fitted Values",
     ylab = "Residuals",
     pch = 19, col = "darkgreen")
abline(h = 0, col = "red", lwd = 2)



# Transport 27 against climatology --------------------------------------------

time_range_trans27 <- transport_data_27 %>%
  reframe(time_range = range(time)) %>%
  pull()


clim_points_trans_27 <- tibble(date = seq(time_range_trans27[1], time_range_trans27[2], by = "1 day")) %>%
  mutate(doy = pmin(yday(date), 365)) %>%
  left_join(climatology_trans_27 %>% select(doy, trans_clim), by = "doy")

clim_points_trans_27 <- clim_points_trans_27 %>% 
  mutate(date = as.Date(date))





trans_with_clim_27 <- transport_data_27 %>%
  mutate(doy = yday(time)) %>%
  left_join(climatology_trans_27 %>% select(doy, trans_clim), by = "doy") %>%
  mutate(
    trans_anomaly = transport_Sv - trans_clim,
    anom_label = case_when(
      trans_anomaly > 0 ~ "Anomaly (+)",
      TRUE ~ "Anomaly (-)")
  )


trans_with_clim_27 <- trans_with_clim_27 %>%
  mutate(date = as.Date(time))


trans_anom_27 <- trans_with_clim_27 %>% 
  select(trans_anomaly)

# Plot strength with climatology
p <- ggplot() +
  # Segments for anomalies
  geom_segment(
    data = trans_with_clim_27,
    aes(x = date, xend = date, y = trans_clim, yend = transport_Sv,
        colour = anom_label, linetype = anom_label),
    linewidth = 1.2
  ) +
  # Climatology line
  geom_line(
    data = clim_points_trans_27,
    aes(x = date, y = trans_clim, linetype = "Climatology", colour = "Climatology"),
    linewidth = 1.1
  ) +
  # Points for observed mean_vcur
  geom_point(
    data = trans_with_clim_27,
    aes(x = date, y = transport_Sv, fill = anom_label, colour = anom_label),
    size = 2.5, shape = 21
  ) +
  scale_colour_manual(
    breaks = c("Anomaly (+)", "Anomaly (-)", "Climatology"),
    values = c("Anomaly (+)" = "red", "Anomaly (-)" = "blue", "Climatology" = "black")
  ) +
  scale_fill_manual(
    breaks = c("Anomaly (+)", "Anomaly (-)"),
    values = c("Anomaly (+)" = "red", "Anomaly (-)" = "blue"),
    guide = "none" 
    
  ) +
  scale_linetype_manual(
    breaks = c("Anomaly (+)", "Anomaly (-)", "Climatology"),
    values = c("solid", "solid", "solid")
  ) +
  scale_x_date(date_breaks = "1 year", minor_breaks = NULL, date_labels = "%Y") +
  #coord_cartesian(ylim = c(-0.3, 0.4)) +
  labs(
    x = "Time",
    y = "Transport (Sv, poleward +)",
    colour = NULL,
    fill = NULL,
    linetype = NULL,
    title = "Current transport Climatology (monthly average)"
  )

ggsave(file.path("output", "trans_with_clim_27.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)



# Transport 32 against climatology --------------------------------------------

time_range_trans32 <- transport_data_32 %>%
  reframe(time_range = range(time)) %>%
  pull()


clim_points_trans_32 <- tibble(date = seq(time_range_trans32[1], time_range_trans32[2], by = "1 day")) %>%
  mutate(doy = pmin(yday(date), 365)) %>%
  left_join(climatology_trans_32 %>% select(doy, trans_clim), by = "doy")

clim_points_trans_32 <- clim_points_trans_32 %>% 
  mutate(date = as.Date(date))





trans_with_clim_32 <- transport_data_32 %>%
  mutate(doy = yday(time)) %>%
  left_join(climatology_trans_32 %>% select(doy, trans_clim), by = "doy") %>%
  mutate(
    trans_anomaly = transport_Sv - trans_clim,
    anom_label = case_when(
      trans_anomaly > 0 ~ "Anomaly (+)",
      TRUE ~ "Anomaly (-)")
  )


trans_with_clim_32 <- trans_with_clim_32 %>%
  mutate(date = as.Date(time))


trans_anom_32 <- trans_with_clim_32 %>% 
  select(trans_anomaly)

# Plot strength with climatology
p <- ggplot() +
  # Segments for anomalies
  geom_segment(
    data = trans_with_clim_32,
    aes(x = date, xend = date, y = trans_clim, yend = transport_Sv,
        colour = anom_label, linetype = anom_label),
    linewidth = 1.2
  ) +
  # Climatology line
  geom_line(
    data = clim_points_trans_32,
    aes(x = date, y = trans_clim, linetype = "Climatology", colour = "Climatology"),
    linewidth = 1.1
  ) +
  # Points for observed mean_vcur
  geom_point(
    data = trans_with_clim_32,
    aes(x = date, y = transport_Sv, fill = anom_label, colour = anom_label),
    size = 2.5, shape = 21
  ) +
  scale_colour_manual(
    breaks = c("Anomaly (+)", "Anomaly (-)", "Climatology"),
    values = c("Anomaly (+)" = "red", "Anomaly (-)" = "blue", "Climatology" = "black")
  ) +
  scale_fill_manual(
    breaks = c("Anomaly (+)", "Anomaly (-)"),
    values = c("Anomaly (+)" = "red", "Anomaly (-)" = "blue"),
    guide = "none" 
    
  ) +
  scale_linetype_manual(
    breaks = c("Anomaly (+)", "Anomaly (-)", "Climatology"),
    values = c("solid", "solid", "solid")
  ) +
  scale_x_date(date_breaks = "1 year", minor_breaks = NULL, date_labels = "%Y") +
  #coord_cartesian(ylim = c(-0.3, 0.4)) +
  labs(
    x = "Time",
    y = "Transport (Sv, poleward +)",
    colour = NULL,
    fill = NULL,
    linetype = NULL,
    title = "Current transport Climatology (monthly average)"
  )

ggsave(file.path("output", "trans_with_clim_32.png"),
       plot = p, width = 1200 / 96, height = 600 / 96, dpi = 96,
       device = png)







