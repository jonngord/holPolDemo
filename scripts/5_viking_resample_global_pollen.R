# This script resamples 150 pollen grains from each sample in each pollen 
# record

# Load pollen standardisation function
source("scripts/standardise_pollen_fun_viking.R")

# Load packages
packages <- c(
  "tidyverse"                       
)
for (package in packages) {
  library(package, character.only = TRUE)
}

# Set locations for outputted resampled datasets to be written to
folder_save <- file.path("data_processed", 
                         paste0("resampled_pollen"))

if ( ! dir.exists(folder_save)) {
  dir.create(folder_save, showWarnings = FALSE)
}

# Read in dataset_ids
global_pollen <- 
  readRDS("data_processed/global_pollen.RDS")

# Read in dataset_ids
dataset_ids <- 
  readRDS("data_processed/global_pollen.RDS") %>% 
  dplyr::pull(dataset_id) %>% 
  unique()

# Choose same 50 dataset_ids
set.seed(123)
samples <- sample(1:length(dataset_ids), 50)

# Resample chosen datasets 10 times (full analysis = 1000 resamples)
# 10 resamples of 50 dataset_ids takes ~ 30 minutes

for (i in 1:10) {   # 10 resamples gives some appreciation of uncertainty but: ---->
#  for (i in 1:2) { # RESAMPLE TWICE TO REDUCE OVERALL COMPUTATION TIME, BUT 
                    # UNCERTAINTY ESTIMATES WILL BE LOST AND FINAL PLOTS WILL
                    # BE DISSIMILAR
  
  # Isolate dataset_id
  ID <- dataset_ids[samples]
  
  # Sort filesave
  file_save <- file.path(folder_save, 
                         paste0("pollen_resample_", i, ".RDS"))
  
  
  ###### ORGANISE DATA ######
  
  # Metadata
  sequence_meta <- 
    global_pollen %>% 
    # Line added for demo -- only doing 50 datasets.
    # Full analayis resamples all datasets 1000 times
    dplyr::filter(dataset_id %in% ID) %>% 
    dplyr::select(event:depth_m)
  
  # Pollen data
  sequence_counts <- 
    global_pollen %>% 
    # Line added for demo -- only doing 50 datasets.
    # Full analayis resamples all datasets 1000 times
    dplyr::filter(dataset_id %in% ID) %>% 
    dplyr::select(-c(event:depth_m)) %>% 
    as.data.frame()
  
  # dataset_id
  sequence_datasetid <- 
    sequence_meta %>% 
    dplyr::select(dataset_id) %>% 
    as.character()
  
  
  ###### RESAMPLE 150 grains from each sample ######
  
  
  # Function from package RRatepol 
  # (Mottl et al, 2021: https://rdrr.io/github/HOPE-UIB-BIO/R-Ratepol-package/src/R/fc_standardise_community_data.R)
  
  standardised_pollen <- 
    fc_standardise_community_data(sequence_counts, 
                                  sequence_datasetid, 
                                  150, 
                                  FALSE)
  
  ###### BIND ######
  standardised_pollen <-
    dplyr::bind_cols(
      sequence_meta,
      standardised_pollen
    )
  
  ###### SAVE ######
  
  saveRDS(object = standardised_pollen, file = file_save)
  
}
  
  
  
  
  