# This script computes the diversity (richness, evenness, turnover [as measured by
# Bray-Curtis]) of the resampled pollen datasets


########################################
######           Config           ######
########################################

# Source functions
source("scripts/bray_curtis_config.R")

packages <- c(
  "tidyverse", "vegan", "broom","glmmTMB", "gtools", "zoo"
)

for (package in packages) {
  
  library(package, character.only = TRUE)
}

# Set location for outputted data to be written to 
folder_save <- "data_processed/diversity_results"

if ( ! dir.exists(folder_save)) {
  dir.create(folder_save, showWarnings = FALSE)
}

# Make directory of 10 pollen resamples 
pollen_resamples <- 
  as.list(
    dir("data_processed/pollen_resample_with_age_draw",
        pattern = ".RDS",
        full.names = TRUE))

# Set folder saves
# Set location for outputted data to be written to 
folder_save <- file.path("data_processed/diversity_results")

if ( ! dir.exists(folder_save)) {
  dir.create(folder_save, showWarnings = FALSE)
}

# Set richness/evenness folder save 
richness_evenness_folder_save <- 
  file.path(folder_save, 
            paste0("richness_evenness_results"))

if ( ! dir.exists(richness_evenness_folder_save)) {
  dir.create(richness_evenness_folder_save, showWarnings = FALSE)
}




# Loop over datasets 
resamples <- 1:length(pollen_resamples)

for (resample in resamples) {
  
  # Load single resample
  global_pollen <- readRDS(pollen_resamples[[resample]])
  
  
  ####### Richness ####### 
  
  
  # Isolate meta data
  meta <- global_pollen %>% 
    dplyr::select(1:16)
  
  # Isolate pollen data
  pollen <- global_pollen %>% 
    dplyr::select(-c(1:16))
  
  # Calculate richness
  pollen_richness <- 
    rowSums(pollen != 0)
  
  # Bind to metadata
  global_pollen_richness <- 
    dplyr::bind_cols(
      meta,
      pollen_richness
    ) %>% 
    rename(richness = 17)
  
  
  
  ####### Evenness ####### 
  
  
  
  H <- vegan::diversity(pollen)      # Shannon index 
  S <- vegan::specnumber(pollen)     # Number of species
  J <- H/log(S)                      # Shannon index/log(species count) = Pielou's evenness
  
  global_pollen_evenness <- 
    dplyr::bind_cols(
      meta,
      J
    ) %>% 
    rename(evenness = 17)
  
  
  ################################################################################
  
  # Save richness and evenness results
  
  # Assign to list and save
  richness_evenness_results <- 
    list(
      global_pollen_richness = global_pollen_richness,
      global_pollen_evenness = global_pollen_evenness
    )
  
  # Set filesave
  file_save <- file.path(richness_evenness_folder_save, 
                         paste0("resample_", resample, "_rich_evenness", '.RDS')) 
  
  ### Save 
  saveRDS(object = richness_evenness_results,
          file = file_save)
  
  # Remove objects and gc()
  rm(global_pollen_evenness,
     global_pollen_richness)
  
  gc()
  
  ################################################################################
  
  
  ####### Turnover ####### 
  
  
  ### Split into dataset_ids
  subset_clean <- 
    global_pollen %>% 
    dplyr::group_by(dataset_id) %>% 
    dplyr::group_split()
  
  ### Compute B/C on all dataframes and output B/C vs. time diffs ###
  
  bray_results <- list()
  for (i in 1:length(subset_clean)) {
    
    # print(i)
    
    # Isolate each dataset_id
    pollen_site <- subset_clean[[i]]
    
    # Run Bray-Curtis computation function on that pollen df
    pollen_site <- 
      pollen_site %>% 
      dplyr::distinct(age_draw, .keep_all = TRUE) %>%  
      dplyr::relocate(age_draw, dplyr::everything()) %>% 
      tidyr::drop_na()
    
    # Isolate ages and id
    pollen_site_age_draw <- pollen_site$age_draw
    pollen_site_dataset_id <- pollen_site$dataset_id[1]
    pollen_site_meta <- pollen_site[, c(1:16)]
    region <- pollen_site_meta$continent[1]
    
    # Calculate midpoint of time for plotting successive diffs 
    age_draw_midpoint <- zoo::rollapply(pollen_site_age_draw, 2, mean)
    
    # Remove metadata from pollen df
    pollen_site <- pollen_site %>% 
      dplyr::select(-c(1:16)) 
    
    # Compute Bray-Curtis
    raw_bc <- vegan::vegdist(pollen_site, method = 'bray')
    raw_bc <- as.data.frame(as.matrix(raw_bc))
    
    # add column and row names
    colnames(raw_bc) <- pollen_site_age_draw
    raw_bc <- cbind(pollen_site_age_draw, raw_bc) %>% 
      dplyr::rename(time = 1)
    
    # Work out the distance between all the rows of the data matrix 
    # (https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/dist)
    time_diff_mat <- as.matrix(dist(raw_bc[,"time"],
                                    diag = TRUE, upper = TRUE))
    
    # Only need differences off the diagonal. So, match matrix sizes 
    # by removing the 'time' column first
    time <- raw_bc[,1] 
    bc_diff <- raw_bc[,-1] 
    
    # Add a row to the time_diff_matrix and check that the matrix columns 
    # and rows have the same dimensions, using '=='.
    delta <- (row(time_diff_mat) + 1) == col(time_diff_mat) 
    
    successive_time_diff <- time_diff_mat[delta] # Time interval matrix
    successive_bc_diff <- bc_diff[delta] # Successive bray-curtis values
    
    # Combine differences
    successive_diffs <- 
      # drop the first, because we need successive differences
      data.frame("time" = raw_bc[-1,"time"], 
                 "bc_diff" = successive_bc_diff,
                 "time_diff" = successive_time_diff,
                 "dataset_id" = pollen_site_dataset_id,
                 "region" = region, 
                 "age_draw_midpoint" = age_draw_midpoint)
    
    # Get all pairwise differences under n years for building the model 
    # (n is the largest time interval between samples)
    # The maximum distance between any 2 samples allowed is 3000 years
    lt <- ((time_diff_mat < 3000) * (time_diff_mat != 0)) == 1 
    all_paired_diffs <- data.frame("time_diff" = time_diff_mat[lt],
                                   "bc_diff" = as.matrix(bc_diff)[lt],
                                   "dataset_id" = pollen_site_dataset_id,
                                   "region" = region)
    
    output <- list(successive_diffs = successive_diffs,
                   all_paired_diffs = all_paired_diffs,
                   age_draw = pollen_site_age_draw)
    
    bray_results[[i]] <- output
    
  }
  
  ### Extract time diff x bc value df
  all_paired_diffs = tibble()
  successive_diffs = tibble()
  
  for (l in seq_along(bray_results)) {
    
    site <- bray_results[[l]]
    time_diffs <- site[["all_paired_diffs"]]
    succ_diffs <- site[["successive_diffs"]]
    all_paired_diffs <- rbind(all_paired_diffs, time_diffs)
    successive_diffs <- rbind(successive_diffs, succ_diffs)
  }
  
  # Recale bc_diff to not include 0 and 1 values
  bc_diff.scaled <- transform01(all_paired_diffs$bc_diff)
  all_paired_diffs$bc_diff.scaled <- bc_diff.scaled
  
  
  ### GLMMTMB FIT
  # Quantify expected amount of bray-curtis for a given time interval
  

  # Run model
  bcdiff_timediff_betareg <- 
    glmmTMB(
      bc_diff.scaled ~ log(time_diff),
      data = all_paired_diffs, 
      family = beta_family(link = "logit"),  
      na.action = "na.omit", 
      dispformula = ~ time_diff,
      control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS"))
    )
  
  summary(bcdiff_timediff_betareg)

  ### Predict from model the expected bc per time interval
  
  preds_population <- 
    data.frame(
      time_diff = all_paired_diffs$time_diff
    )
  
  preds_pop <- predict(bcdiff_timediff_betareg, 
                       newdata = preds_population)
  preds_pop_inv <- plogis(preds_pop) 
  preds_population$preds_population_inverse <- preds_pop_inv
  
  
  ### Predict for sequential samples in dataset actual bc - expected bc for the interval
  ### between the samples
  
  # Make df for regional points
  regional_points <- 
    data.frame(
      time_diff = successive_diffs$time_diff
    )
  
  # Predict
  regional_predictions <- predict(bcdiff_timediff_betareg, 
                                  newdata = regional_points)
  regional_predictions <- plogis(regional_predictions) 
  
  # Add predicted bray-curtis and residuals columns to the dataset
  successive_diffs$predicted <- regional_predictions 
  
  # Subtract actual bray-curtis from expected to adjust for time interval between 
  # pollen samples
  successive_diffs$resid <- (successive_diffs$bc_diff) - successive_diffs$predicted 
  
  ### Combine objects into list for save 
  turnover <- 
    list(
      model = bcdiff_timediff_betareg, 
      all_paired_diffs = all_paired_diffs,
      successive_diffs = successive_diffs
    )
  
  # Set bray folder save 
  bray_folder_save <- 
    file.path(folder_save,
              paste0("turnover_results"))
  
  if ( ! dir.exists(bray_folder_save)) {
    dir.create(bray_folder_save, showWarnings = FALSE)
  }
  
  # Set filesave
  file_save <- file.path(bray_folder_save, 
                         paste0("resample_", resample, "_turnover", '.RDS')) 
  
  ### Save 
  saveRDS(object = turnover,
          file = file_save)
  
}





















