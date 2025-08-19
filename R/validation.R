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

