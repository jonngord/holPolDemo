# This script reads the files outputted by the ARIMA modelling scripts and.
# saves them  to a list 

library(tidyverse)
library(parallel)
library(foreach)
library(doParallel)

# Read arima preds .RDS files into a list
richness <- as.list(dir("data_processed/arima_preds/richness/", 
                        pattern = ".RDS", 
                        full.names = TRUE))

evenness <- as.list(dir("data_processed/arima_preds/evenness/", 
                        pattern = ".RDS", 
                        full.names = TRUE))

turnover <- as.list(dir("data_processed/arima_preds/turnover/", 
                        pattern = ".RDS", 
                        full.names = TRUE))

# Combine into list
all_diversity_preds <- list(
  richness = richness,
  evenness = evenness,
  turnover = turnover
)

# Loop through files and load each one

# Setup parallel backend
cores = detectCores()
clust <- parallel::makeCluster((cores/2) - 1)
registerDoParallel(clust)

preds <- list()

# Load files and extract data
for (diversity_metric in names(all_diversity_preds)) {
  
  files <- all_diversity_preds[[diversity_metric]]
  
  global_arima_preds <-
    
    foreach(i = 1:length(files)) %dopar% {
      
      # Read in each file
      viking_output <- readRDS(files[[i]])
      
      viking_output <- viking_output[["whole_hol_outputs"]]
      
      arima_preds <- viking_output[["weighted_av_preds"]]
      
      predictor_variables <- viking_output[["predictor_variables"]]
      
      resample <- viking_output[["resample"]]
      
      model <- viking_output[["arima_model"]]
      
      # Combine into list
      viking_files <- list(
        resample = resample,
        arima_preds = arima_preds,
        predictor_variables = predictor_variables,
        model = model
      )
      
      return(viking_files)
      
    }
 
  preds[[diversity_metric]] <- global_arima_preds
   
}

stopCluster(clust)


# Save
saveRDS(preds,
        file = "data_processed/arima_preds/arima_preds.RDS")
















