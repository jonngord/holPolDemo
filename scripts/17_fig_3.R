library(tidyverse)
library(patchwork)
library(doParallel)
library(foreach)

`%dopar%` <- foreach::`%dopar%`


##### Data prep #####


# Read the observed diversity .RDS files into a list
richness_evenness <- as.list(dir("data_processed/diversity_results/richness_evenness_results/",
                                 pattern = ".RDS",
                                 full.names = TRUE))

turnover_files <- as.list(dir("data_processed/diversity_results/turnover_results/",
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
          age_draw <= 1100 ~ 1,
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
      dplyr::relocate(slice, dplyr::everything()) %>% 
      dplyr::filter(slice %in% c(1, 10))
    
    # Quantify mean richness of each site in each 1000 year bin
    richness_slice_av <-
      richness_slice %>%
      dplyr::group_by(dataset_id, slice) %>%
      dplyr::summarise(av_richness = median(richness))
    
    
    # Compare proportional richness changes between slice 10 and slice 1
    
    # Reference slice, i.e. slice to which later slices are compared
    # - Top slice
    top_slice <- 10
    
    # - Comparison slice, i.e. slice comparing to reference slice
    bot_slice <- 1
    
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

    dataset_pc_rich_diff <- 
      pc_rich_diff %>% 
      dplyr::mutate(
        pc_rich_diff_bin = dplyr::case_when(
          percent_rich_diff <= - 105                              ~ -110,
          percent_rich_diff >  - 105   & percent_rich_diff <= -95 ~ -100,
          percent_rich_diff >  - 95    & percent_rich_diff <= -85 ~ - 90,
          percent_rich_diff >  - 85    & percent_rich_diff <= -75 ~ - 80,
          percent_rich_diff >  - 75    & percent_rich_diff <= -65 ~ - 70,
          percent_rich_diff >  - 65    & percent_rich_diff <= -55 ~ - 60,
          percent_rich_diff >  - 55    & percent_rich_diff <= -45 ~ - 50,
          percent_rich_diff >  - 45    & percent_rich_diff <= -35 ~ - 40,
          percent_rich_diff >  - 35    & percent_rich_diff <= -25 ~ - 30,
          percent_rich_diff >  - 25    & percent_rich_diff <= -15 ~ - 20,
          percent_rich_diff >  - 15    & percent_rich_diff <= - 5 ~ - 10,
          percent_rich_diff >    -5    & percent_rich_diff <=   5 ~   0,
          percent_rich_diff >     5    & percent_rich_diff <=  15 ~   10,
          percent_rich_diff >    15    & percent_rich_diff <=  25 ~   20,
          percent_rich_diff >    25    & percent_rich_diff <=  35 ~   30,
          percent_rich_diff >    35    & percent_rich_diff <=  45 ~   40,
          percent_rich_diff >    45    & percent_rich_diff <=  55 ~   50,
          percent_rich_diff >    55    & percent_rich_diff <=  65 ~   60,
          percent_rich_diff >    65    & percent_rich_diff <=  75 ~   70,
          percent_rich_diff >    75    & percent_rich_diff <=  85 ~   80,
          percent_rich_diff >    85    & percent_rich_diff <=  95 ~   90,
          percent_rich_diff >    95    & percent_rich_diff <  105 ~  100,
          percent_rich_diff >=  105                               ~  110
          
        )
      ) %>%
      dplyr::group_by(pc_rich_diff_bin) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::mutate(dataset_percent = n / n_sites * 100,
                    mean_pc_diff = mean_pc_diff,
                    median_pc_diff = median_pc_diff,
                    sd_pc_rich_diff = sd_pc,
                    resample = resample,
                    n_sites = n_sites,
                    top_slice = top_slice,
                    bot_slice = bot_slice)
    
    
    ## Evenness
    
    evenness <- file[["global_pollen_evenness"]]
    
    # Attribute slice name to each 1000 years
    evenness_slice <-
      evenness %>%
      dplyr::mutate(
        slice = dplyr::case_when(
          age_draw <= 1100 ~ 1,
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
      dplyr::relocate(slice, dplyr::everything()) %>% 
      dplyr::filter(slice %in% c(1, 10))
    
    # Quantify median evenness of each site in each 1000 year bin
    evenness_slice_av <-
      evenness_slice %>%
      dplyr::group_by(dataset_id, slice) %>%
      dplyr::summarise(av_evenness = median(evenness))
    
    
    # Compare proportional evenness changes between slice 10 and all slice 1
    
    # Reference slice, i.e. slice to which later slices are compared
    # - Top slice
    top_slice <- 10
    
    # - Comparison slice, i.e. slice comparing to reference slice
    bot_slice <- 1
    
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
    
    dataset_pc_evenness_diff <- 
      pc_evenness_diff %>% 
      dplyr::mutate(
        pc_evenness_diff_bin = dplyr::case_when(
          percent_evenness_diff <= - 105                                  ~ -110,
          percent_evenness_diff >  - 105   & percent_evenness_diff <= -95 ~ -100,
          percent_evenness_diff >  - 95    & percent_evenness_diff <= -85 ~ - 90,
          percent_evenness_diff >  - 85    & percent_evenness_diff <= -75 ~ - 80,
          percent_evenness_diff >  - 75    & percent_evenness_diff <= -65 ~ - 70,
          percent_evenness_diff >  - 65    & percent_evenness_diff <= -55 ~ - 60,
          percent_evenness_diff >  - 55    & percent_evenness_diff <= -45 ~ - 50,
          percent_evenness_diff >  - 45    & percent_evenness_diff <= -35 ~ - 40,
          percent_evenness_diff >  - 35    & percent_evenness_diff <= -25 ~ - 30,
          percent_evenness_diff >  - 25    & percent_evenness_diff <= -15 ~ - 20,
          percent_evenness_diff >  - 15    & percent_evenness_diff <= - 5 ~ - 10,
          percent_evenness_diff >    -5    & percent_evenness_diff <=   5 ~   0,
          percent_evenness_diff >     5    & percent_evenness_diff <=  15 ~   10,
          percent_evenness_diff >    15    & percent_evenness_diff <=  25 ~   20,
          percent_evenness_diff >    25    & percent_evenness_diff <=  35 ~   30,
          percent_evenness_diff >    35    & percent_evenness_diff <=  45 ~   40,
          percent_evenness_diff >    45    & percent_evenness_diff <=  55 ~   50,
          percent_evenness_diff >    55    & percent_evenness_diff <=  65 ~   60,
          percent_evenness_diff >    65    & percent_evenness_diff <=  75 ~   70,
          percent_evenness_diff >    75    & percent_evenness_diff <=  85 ~   80,
          percent_evenness_diff >    85    & percent_evenness_diff <=  95 ~   90,
          percent_evenness_diff >    95    & percent_evenness_diff <  105 ~  100,
          percent_evenness_diff >=  105                                   ~  110
        )
      ) %>%
      dplyr::group_by(pc_evenness_diff_bin) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::mutate(dataset_percent = n / n_sites * 100,
                    mean_pc_diff = mean_pc_diff,
                    median_pc_diff = median_pc_diff,
                    sd_pc_evenness_diff = sd_pc,
                    resample = resample,
                    n_sites = n_sites,
                    top_slice = top_slice,
                    bot_slice = bot_slice)

    
    all_outputs <-
      list(
        richness_output = dataset_pc_rich_diff,
        evenness_output = dataset_pc_evenness_diff
      )
    
    return(all_outputs)
    
  }

stopCluster(clust)



## Extract resampled change estimates


## Richness

# Bind resample results tgoether
rich <- lapply(site_changes, function(x){
  x[["richness_output"]]
}) %>%
  dplyr::bind_rows()

# SUmmarise means and medians across resamples
rich_averages <-
  rich %>%
  dplyr::summarise(median_pc_diff = median(median_pc_diff),
                   mean_pc_diff = mean(mean_pc_diff))

# Summarise data across resamples
rich_data <-
  rich %>%
  dplyr::group_by(pc_rich_diff_bin) %>%
  summarise(median_dataset_percent = median(dataset_percent),
            sd_dataset_percent = sd(dataset_percent))  %>%
  dplyr::arrange(desc(pc_rich_diff_bin)) %>% 
  dplyr::mutate(median_pc_diff = rich_averages$median_pc_diff,
                mean_pc_diff = rich_averages$mean_pc_diff)


# Evenness

# Bind resample results tgoether
evenness <- lapply(site_changes, function(x){
  x[["evenness_output"]]
}) %>%
  dplyr::bind_rows()

# SUmmarise means and medians across resamples
evenness_averages <-
  evenness %>%
  dplyr::summarise(median_pc_diff = median(median_pc_diff),
                   mean_pc_diff = mean(mean_pc_diff))

# Summarise data across resamples
evenness_data <-
  evenness %>%
  dplyr::group_by(pc_evenness_diff_bin) %>%
  summarise(median_dataset_percent = median(dataset_percent),
            sd_dataset_percent = sd(dataset_percent))  %>%
  dplyr::arrange(desc(pc_evenness_diff_bin)) %>% 
  dplyr::mutate(median_pc_diff = evenness_averages$median_pc_diff,
                mean_pc_diff = evenness_averages$mean_pc_diff)



# Plots
size <- 17

# Richness
richness_change_plot <-
  ggplot(rich_data, aes(x = pc_rich_diff_bin, y = median_dataset_percent)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  ylab("Sites (%)") +
  xlab("Richness change (%)") +
  ylim(0, 26) +
  xlim(-80, 120) +
  geom_linerange(aes(x = pc_rich_diff_bin,
                     ymin = median_dataset_percent - sd_dataset_percent,
                     ymax = median_dataset_percent + sd_dataset_percent),
                 colour = "orange",
                 alpha = 0.9,
                 size = 1.3) +
  geom_vline(xintercept = 0, colour = "black", size = 1.5) +
  geom_vline(aes(xintercept = mean_pc_diff), colour = "red", linetype = "dashed", size = 1.5) +
  geom_vline(aes(xintercept = median_pc_diff), colour = "red", linetype = "solid", size = 1.5) +
  theme_bw() +
  theme(text = element_text(size=size)) +
  annotate("text", x = -80, y = 25, label = "A", size = size-10)

# Evenness
evenness_change_plot <-
  ggplot(evenness_data, aes(x = pc_evenness_diff_bin, y = median_dataset_percent)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  ylab("Sites (%)") +
  xlab("Evenness change (%)") +
  ylim(0, 26) +
  xlim(-80, 120) +
  geom_linerange(aes(x = pc_evenness_diff_bin,
                     ymin = median_dataset_percent - sd_dataset_percent,
                     ymax = median_dataset_percent + sd_dataset_percent),
                 colour = "orange",
                 alpha = 0.9,
                 size = 1.3) +
  geom_vline(xintercept = 0, colour = "black", size = 1.5, size = 1.5) +
  geom_vline(aes(xintercept = mean_pc_diff), colour = "red", linetype = "dashed", size = 1.5) +
  geom_vline(aes(xintercept = median_pc_diff), colour = "red", linetype = "solid", size = 1.5) +
  theme_bw() +
  theme(text = element_text(size=size), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  annotate("text", x = -80, y = 25, label = "B", size = size-10)



(all_change_plots <-
    richness_change_plot | evenness_change_plot)

# NOTE #

# This analysis doesnt work that well with such a small number of dataset_ids 
# and replicates. To return Fig. 3 from the manuscript requires most of the 
# datasets included

ggsave("figures/fig_03.jpeg",
       plot = all_change_plots,
       width = 12,
       height = 5)

