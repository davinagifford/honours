### read_sst_data.R
###
### Reads SST data for a series of points along the CPR route.
###
### Created: 2025-02-27
### Author: Wayne A. Rochester

library(dplyr)
library(readr)

sites <- read_csv(file.path("data", "sst_sites.csv"), col_types = "idd")
samples_sst <- read_csv(file.path("data", "site_sst.csv"), col_types = "iiTd")

saveRDS(sites, file.path("var", "sst_sites.rds"))
saveRDS(samples_sst, file.path("var", "site_sst.rds"))
