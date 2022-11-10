# This script processes the precipitation data used as a predictor variable in 
# the Holocene diversity ARIMA models 

# Download pastclim package
# devtools::install_github("EvolEcolGroup/pastclim")

library(pastclim)
library(tidyverse)
library(rnaturalearth)
library(ggthemes)
library(sf)
library(foreach)
library(doParallel)
library(Hmisc)


### Code adapted from `Pastclim Applied Example 2` (Leonardi et al, 2022)
# DOI: https://doi.org/10.1101/2022.05.18.492456 

### Dowload Beyer 2020 data
# Data : https://www.nature.com/articles/s41597-020-0552-1
download_dataset(dataset="Beyer2020", bio_variables = "bio12")


# Repeat subsampling procedure 1000 times to arrive at stable 
# global prec. estimate


# time steps of interest
steps <- seq(-18000, 0, by = 1000)

climate_list <- list()
for (step in steps) {
  
  # extract climate (terrestrial surface not covered by ice - but this is fine
  # given the location of pollen core sites)
  climate <- pastclim::region_slice(step, 
                                    "bio12", # precipitation
                                    dataset = "Beyer2020")
  
  climate <- pastclim::df_from_region_slice(climate)
  
  # add time bp
  climate$time_bp <- step
  
  # write into output
  climate_list[[as.character(step)]] <- climate
  
}

# combine them into a single matrix
background_climate <- do.call(rbind, climate_list)
background_climate <- 
  background_climate %>% 
  tibble::rownames_to_column("sample_id")

# Count observations per slice
background_climate %>% 
  dplyr::count(time_bp) 

# Choose only cells that are represented in all time bins
# to get a continuous account
background_climate_filtered <- 
  background_climate %>% 
  dplyr::group_by(x, y) %>% 
  dplyr::mutate(n = n()) %>% 
  dplyr::filter(n == 19)

# Summarise mean prec for each bin
prec_mean <- 
  background_climate_filtered %>% 
  dplyr::group_by(time_bp) %>% 
  dplyr::summarise(mean_prec = mean(bio12, na.rm = TRUE)) %>% 
  tidyr::drop_na()

prec_mean_plot <- 
  ggplot(prec_mean %>% 
           filter(time_bp >= -12000), 
         aes(time_bp, mean_prec)) + 
  geom_line() 

prec <- prec_mean

# Interpolate 
prec <- prec %>% 
  dplyr::mutate(time_bp = seq(18500, 500, -1000))

global_interpol <- 
  tibble::tibble(time_bp = seq(100, 18000, 200)) %>% 
  dplyr::full_join(prec) %>% 
  dplyr::arrange(time_bp)

recent <- 
  global_interpol %>% 
  dplyr::filter(time_bp <= 500) %>% 
  dplyr::mutate(mean_prec = NA)

global_interpol <-
  global_interpol %>% 
  dplyr::filter(time_bp >=  500)

global_interpol$interpol_prec <- zoo::na.approx(global_interpol$mean_prec, 
                                                x = global_interpol$time_bp)  

n <- nrow(global_interpol)

# Calculate change value for precipitation
keep <- tibble(time_bp = seq(500, 18000, 200))
prec_diffs = tibble()
for (i in n:3) {
  
  interpol <- 
    global_interpol %>% 
    dplyr::filter(time_bp %in% keep$time_bp) %>% 
    dplyr::select(time_bp, interpol_prec)
  
  older <- interpol[i, "interpol_prec"]
  younger <- interpol[i-2, "interpol_prec"]
  age <- interpol[i-1, "time_bp"]
  precRoC <- abs(older - younger) # absolute value of the difference between 
  # adjacent bins
  prec_data <- cbind(age, precRoC)
  names(prec_data)[2] <- "prec_change"
  prec_diffs <- rbind(prec_diffs, prec_data)
  
}


# Join prec change and raw prec data
prec_diffs <- 
  prec_diffs %>% 
  dplyr::ungroup() %>% 
  dplyr::bind_rows(recent) %>% 
  dplyr::arrange(time_bp) %>% 
  dplyr::select( - mean_prec) %>% 
  dplyr::full_join(global_interpol,
                   by = "time_bp") %>% 
  dplyr::select(- mean_prec) %>% 
  tidyr::drop_na(prec_change)

# Assign prec change value for final 3 bins to be the same as previous
recent <- tibble(time_bp = c(100, 300, 500), 
                 prec_change = c(prec_diffs[1, 2], prec_diffs[1, 2], prec_diffs[1, 2]),
                 interpol_prec = c(NA, NA, NA))

prec_diffs <- recent %>% dplyr::bind_rows(prec_diffs)

# Plot prec change
mean_change <- 
  ggplot(prec_diffs %>% 
           filter(time_bp <= 12000), aes(time_bp, prec_change)) +
  geom_line() +
  theme_bw() +
  scale_x_reverse(breaks = seq(100, 12000, 500))

# Add leads 
prec_diffs_leads <- 
  prec_diffs %>% 
  dplyr::arrange(time_bp) %>% 
  dplyr::mutate(prec_lead01 = dplyr::lead(interpol_prec, n = 1)) %>%
  dplyr::mutate(prec_lead02 = dplyr::lead(interpol_prec, n = 2)) %>%
  dplyr::mutate(prec_lead03 = dplyr::lead(interpol_prec, n = 3)) %>%
  dplyr::mutate(prec_lead04 = dplyr::lead(interpol_prec, n = 4)) %>%
  dplyr::mutate(prec_lead05 = dplyr::lead(interpol_prec, n = 5)) %>%  
  dplyr::mutate(prec_lead06 = dplyr::lead(interpol_prec, n = 6)) %>%
  dplyr::mutate(prec_lead07 = dplyr::lead(interpol_prec, n = 7)) %>%
  dplyr::mutate(prec_lead08 = dplyr::lead(interpol_prec, n = 8)) %>%
  dplyr::mutate(prec_lead09 = dplyr::lead(interpol_prec, n = 9)) %>%
  dplyr::mutate(prec_lead10 = dplyr::lead(interpol_prec, n = 10)) %>%
  dplyr::mutate(prec_change_lead01 = dplyr::lead(prec_change , n = 1)) %>% 
  dplyr::mutate(precd_change_lead02 = dplyr::lead(prec_change , n = 2)) %>% 
  dplyr::mutate(prec_change_lead03 = dplyr::lead(prec_change , n = 3)) %>% 
  dplyr::mutate(prec_change_lead04 = dplyr::lead(prec_change , n = 4)) %>% 
  dplyr::mutate(prec_change_lead05 = dplyr::lead(prec_change , n = 5)) %>% 
  dplyr::mutate(prec_change_lead06 = dplyr::lead(prec_change , n = 6)) %>% 
  dplyr::mutate(prec_change_lead07 = dplyr::lead(prec_change , n = 7)) %>% 
  dplyr::mutate(prec_change_lead08 = dplyr::lead(prec_change , n = 8)) %>% 
  dplyr::mutate(prec_change_lead09 = dplyr::lead(prec_change , n = 9)) %>% 
  dplyr::mutate(prec_change_lead10 = dplyr::lead(prec_change , n = 10)) %>% 
  dplyr::relocate(time_bp, everything()) %>% 
  dplyr::rename(prec = interpol_prec) %>% 
  dplyr::relocate(prec_change, .after = prec_lead10)


saveRDS(prec_diffs_leads, "data_processed/global_prec_leads.RDS")




