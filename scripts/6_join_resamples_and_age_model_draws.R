# This script joins a draw from the posterior distributions of the 
# Bchron age-depth model with one of the 1:1000 resampled global pollen datasets
# Script run on the University of York's researcher servers

# Load libs etc
source("scripts/load_libraries.R")

`%dopar%` <- foreach::`%dopar%`

# Read in top/bottom slice depths for part seqs
slice_depths <- 
  readRDS("data_processed/top_bottom_depths_part_sequences.RDS")

# Read in global_pollen.RDS
global_pollen <- readRDS("data_processed/global_pollen.RDS")

# 1 - Make directory of 10 pollen resamples 
pollen_files <- 
  as.list(
    dir("data_processed/resampled_pollen",
        pattern = ".RDS",
        full.names = TRUE))


# 2 - Read in age-depth models 
age_depth_model_files <- 
  as.list(
    dir("data_processed/age_depths_bchron",
        pattern = ".RDS",
        full.names = TRUE))

age_depth_models = list()
for (f in seq_along(age_depth_model_files)) {
  
  age_model <- readRDS(age_depth_model_files[[f]])
  age_model <- age_model[["bchron_age_depth"]]
  age_depth_models[[f]] <- age_model
  
}

age_depth_models <- 
  age_depth_models %>% 
  dplyr::bind_rows()

# Set location for outputted data to be written to 
folder_save <- "data_processed/pollen_resample_with_age_draw"

if ( ! dir.exists(folder_save)) {
  dir.create(folder_save, showWarnings = FALSE)
}

# 3 - Read in and join 1:10 age-depth draws to 
#     each of the 1:10 pollen datasets

# Setup parallel backend
clust <- parallel::makeCluster(4)
doParallel::registerDoParallel(clust)

out <- 
  
  foreach::foreach(i =  1:10) %dopar%  {
    
    library(tidyverse)
    
    # Read in each pollen resample
    # join to metadata
    # convert depths (currently in m) to cm
    global_pollen_resample <- 
      readRDS(pollen_files[[i]]) %>%              
      dplyr::relocate(data_source:depth_m, dplyr::everything()) %>% 
      dplyr::mutate(depths = depth_m * 100) %>%                       # convert to cm
      dplyr::relocate(depths, .after = depth_m) %>% 
      dplyr::add_count(dataset_id, depths) %>%                        
      dplyr::filter(!n > 1) %>%                                       
      dplyr::select(- n)
    
    # Convert columns to characters (using paste) for join below
    global_pollen_resample$depths <- paste(global_pollen_resample$depths)
    global_pollen_resample$dataset_id <- paste(global_pollen_resample$dataset_id)
    
    # Isolate age model dataset_id and depth
    global_age_depth_model_draw_meta <- 
      age_depth_models %>% 
      dplyr::select(dataset_id, depths)
    
    # Isolate age model ages
    global_age_depth_model_draw_ages <- 
      age_depth_models %>% 
      dplyr::select(- c(dataset_id, depths, median_age) )
    
    # Isolate single draw from age-depth model
    age_draw <- global_age_depth_model_draw_ages[, i]
    global_age_depth_model_draw_meta$age_draw <- age_draw
    
    # Convert columns to characters (using paste) for join below
    global_age_depth_model_draw_meta$dataset_id <- paste(global_age_depth_model_draw_meta$dataset_id)
    global_age_depth_model_draw_meta$depths <- paste(global_age_depth_model_draw_meta$depths)
    
    # Join age draw and pollen resample
    pollen_age_draw <-
      dplyr::left_join(
        global_age_depth_model_draw_meta,
        global_pollen_resample,
        by = c(
          "dataset_id" = "dataset_id",
          "depths" = "depths"
        ) 
      ) %>% 
      dplyr::relocate(age_draw, .before = depth_m) %>% 
      dplyr::mutate(depths = as.numeric(depths)) 
    
    # Slice top and bottom depths of part pollen seqs
    # Isolate full and part seqs
    # Full seqs
    full_seqs <- 
      pollen_age_draw %>% 
      dplyr::filter(!dataset_id %in% slice_depths$dataset_id)
    
    # Part seqs
    part_seqs <-
      pollen_age_draw %>% 
      dplyr::filter(dataset_id %in% slice_depths$dataset_id)
    
    # Slice part seqs by depths given in slice_depths df
    part_sequences_sliced <- list()
    for(ID in unique(part_seqs$dataset_id)){
      
      dataset_id_depth <- 
        slice_depths %>% 
        dplyr::filter(dataset_id == ID)
      
      top_depth <- dataset_id_depth$topdepth
      bottom_depth <- dataset_id_depth$botdepth
      
      dataset_id_part_seq <- 
        part_seqs %>% 
        dplyr::filter(
          dataset_id == ID,
          depths > top_depth & depths < bottom_depth
        )
      
      part_sequences_sliced[[ID]] <- dataset_id_part_seq
      
    }
    
    part_sequences_sliced <- 
      part_sequences_sliced %>% 
      dplyr::bind_rows()
    
    # Bind with full_seqs to recombine global dataset
    pollen_full_part_joined <-
      dplyr::bind_rows(
        full_seqs,
        part_sequences_sliced
      )
    
    # Filter out: 
    # 1 - pre-Holocene samples 
    # 2 - post-1850 samples
    global_pollen_final <- 
      pollen_full_part_joined %>% 
      dplyr::filter(age_draw >= 100 &
                      age_draw <= 11700)
    
    # ID datasets with < 5 samples
    sample_n <-
      global_pollen_final %>% 
      dplyr::group_by(dataset_id) %>% 
      dplyr::summarise(samples = n()) %>% 
      dplyr::filter(samples < 5)
    
    # Filter global dataset to NOT include these samples ^
    global_pollen_final <-
      global_pollen_final %>% 
      dplyr::filter(!dataset_id %in% sample_n$dataset_id)
    
    # save
    filesave <- paste0(folder_save, "/resample_", i, ".RDS")
    
    saveRDS(
      object = global_pollen_final,
      file = filesave
    )
    
    return(i)
    
  }









