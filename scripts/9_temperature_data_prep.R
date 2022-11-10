# This script processes the temperature data used as a predictor variable in 
# the Holocene diversity ARIMA models                

library(tidyverse)


###### Palaeoclimate data ######

# Last 24000 years of ST data from Osman et al (2021) 

# Download from ------>
# DOI: https://doi.org/10.1038/s41586-021-03984-4
# and put in `data_raw`

# Load data, rename columns
lgmr_st <- read_csv("data_raw/osman_etal_2021_LGMR.csv")
lgmr_st <- lgmr_st %>% 
  dplyr::rename(median = `50th`) %>% 
  dplyr::mutate(bin_age_centre = seq(100, 23900, 200)) %>% 
  dplyr::relocate(bin_age_centre, everything())


# Calculate the mean rate of climate change over the 200 year period, 
# in relation to the amount of climate change over a ~400 year period 
# (centred on each 200 yr bin; abs difference between previous and subsequent 
# climate bins) 
climate_diffs = tibble()
for (i in 120:3) {
  
  older <- lgmr_st[i, "median"]
  younger <- lgmr_st[i-2, "median"]
  age <- lgmr_st[i-1, "bin_age_centre"]
  climateRoC <- abs(older - younger) # absolute value of the difference between 
  # adjacent bins
  climate <- cbind(age, climateRoC)
  names(climate)[2] <- "gmst_change"
  climate_diffs <- rbind(climate_diffs, climate)
  
}

# Add climate leads
# Add leads (opposite of lag) of different durations to Holocene data 
climate_leads <- climate_diffs %>% 
  dplyr::arrange(bin_age_centre) %>% 
  dplyr::mutate(gmst_change_lead1 = dplyr::lead(gmst_change, n = 1)) %>% 
  dplyr::mutate(gmst_change_lead2 = dplyr::lead(gmst_change, n = 2)) %>% 
  dplyr::mutate(gmst_change_lead3 = dplyr::lead(gmst_change, n = 3)) %>% 
  dplyr::mutate(gmst_change_lead4 = dplyr::lead(gmst_change, n = 4)) %>%
  dplyr::mutate(gmst_change_lead5 = dplyr::lead(gmst_change, n = 5)) %>%
  dplyr::mutate(gmst_change_lead6 = dplyr::lead(gmst_change, n = 6 )) %>%
  dplyr::mutate(gmst_change_lead7 = dplyr::lead(gmst_change, n = 7)) %>%
  dplyr::mutate(gmst_change_lead8 = dplyr::lead(gmst_change, n = 8 )) %>%
  dplyr::mutate(gmst_change_lead9 = dplyr::lead(gmst_change, n = 9)) %>%
  dplyr::mutate(gmst_change_lead10 = dplyr::lead(gmst_change, n = 10)) %>%
  dplyr::mutate(gmst_change_lead11 = dplyr::lead(gmst_change, n = 11)) %>%
  dplyr::mutate(gmst_change_lead12 = dplyr::lead(gmst_change, n = 12)) %>%
  dplyr::mutate(gmst_change_lead13 = dplyr::lead(gmst_change, n = 13)) %>%
  dplyr::mutate(gmst_change_lead14 = dplyr::lead(gmst_change, n = 14)) %>%
  dplyr::mutate(gmst_change_lead15 = dplyr::lead(gmst_change, n = 15)) 


saveRDS(climate_leads, "data_processed/global_binned_climate_data.RDS")
