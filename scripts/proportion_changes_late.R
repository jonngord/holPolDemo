library(tidyverse)
library(patchwork)
library(doParallel)
library(foreach)
library(gt)

`%dopar%` <- foreach::`%dopar%`


##### Data prep #####


# Read the observed diversity .RDS files into a list
richness_evenness <- as.list(dir("data_processed/diversity_results/richness_evenness_results/", 
                                 pattern = ".RDS", 
                                 full.names = TRUE))


# Loop through files and load each one

# Setup parallel backend
cores = detectCores()
clust <- parallel::makeCluster((cores/2))
registerDoParallel(clust)

site_changes <-
  
  foreach(ii = 1:length(richness_evenness)) %dopar% {
    
    library(tidyverse)
    
    # Set array index
    resample <- ii
    
    # Read in richness and evenness files 
    file <- readRDS(richness_evenness[[resample]])
    
    ##### Change analysis #####
    
    ## Richness
    
    richness <- file[["global_pollen_richness"]]
    
    # Attribute slice name to each 1000 years
    richness_slice <- 
      richness %>% 
      dplyr::mutate(
        slice = dplyr::case_when(
          age_draw <= 1000 ~ 1,
          age_draw > 1000 & age_draw <= 2000 ~ 2,
          age_draw > 2000 & age_draw <= 3000 ~ 3,
          age_draw > 3000 & age_draw <= 4000 ~ 4,
          age_draw > 4000 & age_draw <= 5000 ~ 5,
          age_draw > 5000 & age_draw <= 6000 ~ 6,
          age_draw > 6000 & age_draw <= 7000 ~ 7,
          age_draw > 7000 & age_draw <= 8000 ~ 8,
          age_draw > 8000 & age_draw <= 9000 ~ 9,
          age_draw > 9000 & age_draw <= 10000 ~ 10, 
          age_draw > 10000 & age_draw <= 11000 ~ 11, 
          age_draw > 11000 & age_draw <= 12000 ~ 12
        )
      ) %>% 
      dplyr::relocate(slice, dplyr::everything())
    
    # Quantify mean richness of each site in each 1000 year bin
    richness_slice_av <- 
      richness_slice %>% 
      dplyr::group_by(dataset_id, slice) %>% 
      dplyr::summarise(av_richness = median(richness))
    
    
    # Compare proportional richness changes between slice 10 and all subsequent 
    # slices
    richness_post_10k <- data.frame()
    for(i in 9:1) {
      
      # Reference slice, i.e. slice to which later slices are compared
      # - Top slice
      top_slice <- 10
      
      # - Comparison slice, i.e. slice comparing to reference slice
      bot_slice <- i
      
      # Richness values from top slice (lowest richness slice = 9000 cal. years BP
      # and therefore is our comparison)
      top_slice_richness <- 
        richness_slice_av %>% 
        dplyr::filter(slice == top_slice) %>% 
        dplyr::rename(av_richness_reference_slice = 3, 
                      reference_slice = slice)
      
      # dataset_ids contained in this slice
      top_slice_ids <- 
        top_slice_richness %>% 
        dplyr::pull(dataset_id) %>% 
        unique()
      
      # Richness values from bottom slice from dataset_ids that were also present in
      # the top (earliest) slice
      bot_slice_richness <-
        richness_slice_av %>% 
        dplyr::filter(slice == bot_slice, 
                      dataset_id %in% top_slice_ids) %>% 
        dplyr::rename(av_richness_comparison_slice = 3,
                      comparison_slice = slice)
      
      
      # Right join to select only sites that were present in the bottom slice as well 
      # as the top 
      top_bot_richness_comp <- 
        dplyr::right_join(
          top_slice_richness, 
          bot_slice_richness, 
          by = "dataset_id"
        )
      
      n_sites <- length(unique(top_bot_richness_comp$dataset_id))
      
      # Calculate and plot proportional change
      pc_rich_diff <- 
        top_bot_richness_comp %>% 
        dplyr::mutate(rich_diff = av_richness_comparison_slice - av_richness_reference_slice,
                      percent_rich_diff = rich_diff / av_richness_reference_slice * 100) %>% 
        tidyr::drop_na(percent_rich_diff) 
      
      n_sites <- nrow(pc_rich_diff)
      mean_pc_diff <- mean(pc_rich_diff$percent_rich_diff)
      median_pc_diff <- median(pc_rich_diff$percent_rich_diff)
      sd_pc <- sd(pc_rich_diff$percent_rich_diff)
      
      total_sites_w_positive_changes <- 
        pc_rich_diff %>% 
        dplyr::filter(rich_diff > 0) %>% 
        nrow()
      
      total_sites_w_positive_changes_per_cent <- total_sites_w_positive_changes / n_sites * 100
      
      dataset_pc_rich_diff <-  
        pc_rich_diff %>% 
        dplyr::mutate(pc_rich_diff_bin = floor(percent_rich_diff / 10) * 10) %>% 
        dplyr::group_by(pc_rich_diff_bin) %>% 
        dplyr::summarise(n = n()) %>% 
        dplyr::mutate(dataset_percent = n / n_sites * 100,
                      mean_pc_diff = mean_pc_diff,
                      median_pc_diff = median_pc_diff,
                      sd_pc_rich_diff = sd_pc, 
                      resample = resample,
                      n_sites = n_sites, 
                      top_slice = top_slice,
                      bot_slice = bot_slice, 
                      pc_positive = total_sites_w_positive_changes_per_cent)
      
      
      
      richness_post_10k <- dplyr::bind_rows(richness_post_10k, dataset_pc_rich_diff)
      
    }
    
    ## Evenness
    
    
    evenness <- file[["global_pollen_evenness"]]
    
    # Attribute slice name to each 1000 years
    evenness_slice <- 
      evenness %>% 
      dplyr::mutate(
        slice = dplyr::case_when(
          age_draw <= 1000 ~ 1,
          age_draw > 1000 & age_draw <= 2000 ~ 2,
          age_draw > 2000 & age_draw <= 3000 ~ 3,
          age_draw > 3000 & age_draw <= 4000 ~ 4,
          age_draw > 4000 & age_draw <= 5000 ~ 5,
          age_draw > 5000 & age_draw <= 6000 ~ 6,
          age_draw > 6000 & age_draw <= 7000 ~ 7,
          age_draw > 7000 & age_draw <= 8000 ~ 8,
          age_draw > 8000 & age_draw <= 9000 ~ 9,
          age_draw > 9000 & age_draw <= 10000 ~ 10, 
          age_draw > 10000 & age_draw <= 11000 ~ 11, 
          age_draw > 11000 & age_draw <= 12000 ~ 12
        )
      ) %>% 
      dplyr::relocate(slice, dplyr::everything())
    
    # Quantify median evenness of each site in each 1000 year bin
    evenness_slice_av <- 
      evenness_slice %>% 
      dplyr::group_by(dataset_id, slice) %>% 
      dplyr::summarise(av_evenness = median(evenness))
    
    
    # Compare proportional evenness changes between slice 10 and all subsequent 
    # slices
    evenness_post_10k <- data.frame()
    for(i in 9:1) {
      
      # Reference slice, i.e. slice to which later slices are compared
      # - Top slice
      top_slice <- 10
      
      # - Comparison slice, i.e. slice comparing to reference slice
      bot_slice <- i
      
      # evenness values from top slice (lowest evenness slice = 9000 cal. years BP
      # and therefore is our comparison)
      top_slice_evenness <- 
        evenness_slice_av %>% 
        dplyr::filter(slice == top_slice) %>% 
        dplyr::rename(av_evenness_reference_slice = 3, 
                      reference_slice = slice)
      
      # dataset_ids contained in this slice
      top_slice_ids <- 
        top_slice_evenness %>% 
        dplyr::pull(dataset_id) %>% 
        unique()
      
      # evenness values from bottom slice from dataset_ids that were also present in
      # the top (earliest) slice
      bot_slice_evenness <-
        evenness_slice_av %>% 
        dplyr::filter(slice == bot_slice, 
                      dataset_id %in% top_slice_ids) %>% 
        dplyr::rename(av_evenness_comparison_slice = 3,
                      comparison_slice = slice)
      
      
      # Right join to select only sites that were present in the bottom slice as well 
      # as the top 
      top_bot_evenness_comp <- 
        dplyr::right_join(
          top_slice_evenness, 
          bot_slice_evenness, 
          by = "dataset_id"
        )
      
      n_sites <- length(unique(top_bot_evenness_comp$dataset_id))
      
      # Calculate and plot proportional change
      pc_evenness_diff <- 
        top_bot_evenness_comp %>% 
        dplyr::mutate(evenness_diff = av_evenness_comparison_slice - av_evenness_reference_slice,
                      percent_evenness_diff = evenness_diff / av_evenness_reference_slice * 100) %>% 
        tidyr::drop_na(percent_evenness_diff) 
      
      n_sites <- nrow(pc_evenness_diff)
      mean_pc_diff <- mean(pc_evenness_diff$percent_evenness_diff)
      median_pc_diff <- median(pc_evenness_diff$percent_evenness_diff)
      sd_pc <- sd(pc_evenness_diff$percent_evenness_diff)
      
      total_sites_w_positive_changes <- 
        pc_evenness_diff %>% 
        dplyr::filter(evenness_diff > 0) %>% 
        nrow()
      
      total_sites_w_positive_changes_per_cent <- total_sites_w_positive_changes / n_sites * 100
      
      dataset_pc_evenness_diff <-  
        pc_evenness_diff %>% 
        dplyr::mutate(pc_evenness_diff_bin = floor(percent_evenness_diff / 10) * 10) %>% 
        dplyr::group_by(pc_evenness_diff_bin) %>% 
        dplyr::summarise(n = n()) %>% 
        dplyr::mutate(dataset_percent = n / n_sites * 100,
                      mean_pc_diff = mean_pc_diff,
                      median_pc_diff = median_pc_diff,
                      sd_pc_evenness_diff = sd_pc, 
                      resample = resample,
                      n_sites = n_sites, 
                      top_slice = top_slice,
                      bot_slice = bot_slice, 
                      pc_positive = total_sites_w_positive_changes_per_cent)
      
      evenness_post_10k <- dplyr::bind_rows(evenness_post_10k, dataset_pc_evenness_diff)
      
    }
    
    
    # Save
    
    all_outputs <- 
      list(
        richness_output = richness_post_10k,
        evenness_output = evenness_post_10k
      )
    
    return(all_outputs)
    
  }

stopCluster(clust)


## Extract resampled change estimates 

## Richness

# Bind resample results tgoether
rich_post <- lapply(site_changes, function(x){
  x[["richness_output"]]
}) %>% 
  dplyr::bind_rows() 

# Summarise data across resamples
richness_positive_changes_average <- 
  rich_post %>% 
  dplyr::group_by(bot_slice) %>% 
  summarise(median_positive_pc = round(median(pc_positive)),
            mean_positive_pc = round(mean(pc_positive)),
            sd_dataset_percent = round(sd(pc_positive))) %>% 
  dplyr::arrange(desc(bot_slice)) %>% 
  dplyr::rename(rich_median_positive_pc = median_positive_pc) %>% 
  dplyr::select(bot_slice, rich_median_positive_pc)


# Evenness

# Post (bin 10)

# Bind resample results tgoether
evenness_post <- lapply(site_changes, function(x){
  x[["evenness_output"]]
}) %>% 
  dplyr::bind_rows()

# Summarise data across resamples
evenness_positive_changes_average <- 
  evenness_post %>% 
  dplyr::group_by(bot_slice) %>% 
  summarise(median_positive_pc = round(median(pc_positive)),
            mean_positive_pc = round(mean(pc_positive)),
            sd_dataset_percent = round(sd(pc_positive))) %>% 
  dplyr::rename(evenness_median_positive_pc = median_positive_pc) %>% 
  dplyr::select(bot_slice, evenness_median_positive_pc)

both <- richness_positive_changes_average %>% 
  dplyr::left_join(evenness_positive_changes_average, by = "bot_slice")

# Make tables for supp mats

slices <- c("9000-8000","8000-7000", "7000-6000", "6000-5000", 
            "5000-4000", "4000-3000", "3000-2000", "2000-1000",
            "1000-100")

both %>%
  dplyr::select(bot_slice, rich_median_positive_pc, evenness_median_positive_pc) %>% 
  dplyr::mutate(slice = slices) %>% 
  dplyr::relocate(slice, .after = bot_slice) %>%
  dplyr::arrange(bot_slice)  %>% 
  gt::gt() %>% 
  gt::cols_label(
    slice = md("Age range (yBP)"),
    bot_slice = md("Comparison time bin"),
    rich_median_positive_pc = md("Sites with positive richness change (%)"), 
    evenness_median_positive_pc = md("Sites with positive evenness change (%)")) %>%
  cols_align(
    align = "center"
  ) %>% 
gt::gtsave("figures/prop_changes.png", expand = 2)

# NOTE #

# Table will differ significantly from manuscript version with low number of dataset_ids
# and replicates included in analysis -- this is just to demonstrate method


