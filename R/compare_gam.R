# Compare performance of differnt GAMs
# 
# 

# 19 May 2025
# DAvina Gifford

# extract eac_cci value for original GAM

eac_cci_orig <- samples %>% 
  select(eac_cc_original = eac_cci) 

# extract eac_cci value for new GAM

eac_cci_new <- samples %>% 
  select(eac_cci_new = eac_cci)

# compare the two, to see if there is much difference
cor(as.numeric(unlist(eac_cci_orig)), as.numeric(unlist(eac_cci_new)), use = "complete.obs")
 