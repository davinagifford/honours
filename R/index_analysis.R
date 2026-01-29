## index_analysis.R
##
## Analysis of results 
##
## creation date:
## created by: Davina Gifford
## last updated: 2026-01-29
## updated by: Davina Gifford

## --------------------------------------------------------------------------------
# load libraries  
library(tidyverse) 
library(marginaleffects) 
library(gglm)


## Raw: EAC strength vs. copepod index--------------------------------------------------------------------------------


# data wrangling

# Extract data sources
raw_data <- readRDS(file.path("var", "raw_data_list.rds"))

raw_strength_data <- raw_data$raw_strength_data
raw_index_data <- raw_data$raw_index_data
str_clim <- raw_data$str_clim

# combine raw data values

comb_raw_data <- raw_strength_data %>% 
  full_join(raw_index_data, by = "date") %>%  
  mutate(
    month = month(date) # extract month value
  ) %>% 
  rename(mean_strength = mean_vel) %>% 
  filter(!is.na(mean_strength) & !is.na(eac_cci))

comb_raw_data <- comb_raw_data %>% 
  inner_join(str_clim, by = "month")


head(comb_raw_data)




## --------------------------------------------------------------------------------

# create the model

cci_str_model <- lm(eac_cci ~ mean_strength, data = comb_raw_data)

summary(cci_str_model)



## --------------------------------------------------------------------------------
# Index against strength climatology
cci_str_model2 <- lm(eac_cci ~ month_clim_str, data = comb_raw_data)
# month_clim_str is the strength climatology value.

summary(cci_str_model2)



## --------------------------------------------------------------------------------
# Index against strength climatology
cci_str_model3 <- lm(eac_cci ~ month_clim_str + mean_strength, data = comb_raw_data)

summary(cci_str_model3)


cci_str_model4 <- lm(mean_strength ~ eac_cci, data = comb_raw_data)

summary(cci_str_model4)


## Model comparison --------------------------------------------------------------------------------



# We test the significance of the non-seasonal component with ANOVA of the nested models.

# compare the two models 

anova(cci_str_model2, cci_str_model3, test = "F")


## Model diagnostics --------------------------------------------------------------------------------

# look at model 3

gglm(cci_str_model3)



plot_predictions(cci_str_model3, condition = "mean_strength") +
  geom_point(data = comb_raw_data,
             aes(x = mean_strength, y = eac_cci))

plot_predictions(cci_str_model3, condition = "month_clim_str") +
  geom_point(data = comb_raw_data,
             aes(x = month_clim_str, y = eac_cci))

plot_predictions(cci_str_model3, condition = "mean_strength", points= 0.25)

plot_predictions(cci_str_model3, condition = "month_clim_str", points= 0.25)


## Mooring Array --------------------------------------------------------------------------------
# load mooring array data

mooring_data <- readRDS(file.path("var", "mooring_data_list.rds"))

mooring_vel_data <- mooring_data$mooring_vel_data 
moor_clim <- mooring_data$moor_clim

comb_moor_data <- mooring_vel_data %>% 
   full_join(moor_clim, by = "date") %>% 
  mutate(
    month = month(date)
  )

head(comb_moor_data)


#A linear regression was conducted to examine the relationship between the Index values and the observed mean monthly velocity.

# create the model

cci_moor_model <- lm(eac_cci ~ mean_vel, data = comb_moor_data)

summary(cci_moor_model)

# A second regression was conducted to understand the seasonal variation.

# Index against velocity climatology 
cci_moor_model2 <- lm(eac_cci ~ vel_clim, data = comb_moor_data) # vel_clim is the  climatology value.  

summary(cci_moor_model2) 


# total variation in CCI that is explained by the seasonal and total variation in curren
# Index against velocity and climatology 
cci_moor_model3 <- lm(eac_cci ~ vel_clim + mean_vel, data = comb_moor_data)  

summary(cci_moor_model3)


## Observed Velocity against the EAC CCI from the Northern part of the EAC --------------------------------------------------------------------------------
# load northern Index values 
cci_north_data <- readRDS(file.path("var", "eac_cci_north_list.rds"))

cci_north <- cci_north_data$cci_north

comb_moor_data_north <- comb_moor_data %>%  
  full_join(cci_north, by = "date")

head(comb_moor_data_north)



# initial model 

cci_moor_model_n <- lm(cci_north ~ mean_vel, data = comb_moor_data_north)

summary(cci_moor_model_n)


# understand the seasonal variation
# Index against velocity climatology  
cci_moor_model_n2 <- lm(cci_north ~ vel_clim, data = comb_moor_data_north) # vel_clim is the  climatology value.    

summary(cci_moor_model_n2) 


# total variation in CCI that is explained by the seasonal and total variation in current
# Index against velocity and climatology  
cci_moor_model_n3 <- lm(cci_north ~ vel_clim + mean_vel, data = comb_moor_data_north)    
summary(cci_moor_model_n3)


## Model comparison --------------------------------------------------------------------------------
# compare the two models   
anova(cci_moor_model_n2, cci_moor_model_n3, test = "F")


## Model diagnostics --------------------------------------------------------------------------------
# look at model 1  
gglm(cci_moor_model_n2)  



plot_predictions(cci_moor_model_n2, condition = "vel_clim") +   geom_point(data = comb_moor_data_north,                                                aes(x = vel_clim, y = cci_north))


## North vs South --------------------------------------------------------------------------------
# load north and south data
cci_data_n <- readRDS(file.path("var", "eac_cci_north.rds"))


index_n <- cci_data_n$month_data %>% # monthly average Index value
  mutate(date = as.Date(trip_month))

cci_data_s <- readRDS(file.path("var", "eac_cci_south.rds"))


index_s <- cci_data_s$month_data %>% # monthly average Index value
  mutate(date = as.Date(trip_month))

# Extract Velocity data for North and South

# Extract data sources
raw_data_n <- readRDS(file.path("var", "raw_data_list_n.rds"))

raw_strength_data_n <- raw_data_n$raw_strength_data_n
str_clim_n <- raw_data_n$str_clim_n

# combine raw data values

comb_raw_data_n <- raw_strength_data_n %>% 
  full_join(index_n, by = "date") %>%  
  mutate(
    month = month(date) # extract month value
  ) %>% 
  rename(mean_strength = mean_vel) %>% 
  filter(!is.na(mean_strength) & !is.na(eac_cci))

comb_raw_data_n <- comb_raw_data_n %>% 
  inner_join(str_clim_n, by = "month")


head(comb_raw_data_n)

# Extract data sources
raw_data_s <- readRDS(file.path("var", "raw_data_list_s.rds"))

raw_strength_data_s <- raw_data_s$raw_strength_data_s
str_clim_s <- raw_data_s$str_clim_s

# combine raw data values

comb_raw_data_s <- raw_strength_data_s %>% 
  full_join(index_s, by = "date") %>%  
  mutate(
    month = month(date) # extract month value
  ) %>% 
  rename(mean_strength = mean_vel) %>% 
  filter(!is.na(mean_strength) & !is.na(eac_cci))

comb_raw_data_s <- comb_raw_data_s %>% 
  inner_join(str_clim_s, by = "month")


head(comb_raw_data_s)





## North of the Separation Zone--------------------------------------------------------------------------------
# create the model

cci_str_model_n <- lm(eac_cci ~ mean_strength, data = comb_raw_data_n)

summary(cci_str_model_n)


# understand the seasonal variation

# Index against velocity climatology 
cci_str_model_n2 <- lm(eac_cci ~ month_clim_str, data = comb_raw_data_n)

summary(cci_str_model_n2)


# total variation in CCI that is explained by the seasonal and total variation in current
cci_str_model_n3 <- lm(eac_cci ~ month_clim_str + mean_strength, data = comb_raw_data_n)

summary(cci_str_model_n3)


## Model comparison --------------------------------------------------------------------------------
# compare the models 

anova(cci_str_model_n2, cci_str_model_n3, test = "F")


## Model diagnostics --------------------------------------------------------------------------------
# look at model 3  
gglm(cci_str_model_n3)  



plot_predictions(cci_str_model_n3, condition = "mean_strength") +
  geom_point(data = comb_raw_data_n,
             aes(x = mean_strength, y = eac_cci))

plot_predictions(cci_str_model_n3, condition = "mean_strength", points= 0.25)

plot_predictions(cci_str_model_n3, condition = "month_clim_str", points= 0.25)


## South of the separation zone--------------------------------------------------------------------------------
# create the initial model  
cci_str_model_s <- lm(eac_cci ~ mean_strength, data = comb_raw_data_s)  
summary(cci_str_model_s)


# understand the seasonal variation
# Index against velocity climatology  
cci_str_model_s2 <- lm(eac_cci ~ month_clim_str, data = comb_raw_data_s)  
summary(cci_str_model_s2)


# total variation in CCI that is explained by the seasonal and total variation in current
cci_str_model_s3 <- lm(eac_cci ~ month_clim_str + mean_strength, data = comb_raw_data_s)  
summary(cci_str_model_s3)




## Model Comparison --------------------------------------------------------------------------------
# compare the models   
anova(cci_str_model_s2, cci_str_model_s3, test = "F")


## Model diagnostics --------------------------------------------------------------------------------
# look at model 3   
gglm(cci_str_model_s3)  


plot_predictions(cci_str_model_s3, condition = "mean_strength") +   
  geom_point(data = comb_raw_data_s,
             aes(x = mean_strength, y = eac_cci))

