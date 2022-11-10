# This script identifies the top and bottom depths at which to slice records 
# that only partially met the inclusion criteria (based on chron. sampling)


source("scripts/load_libraries")

# Read in pollen chronology .RDS
chrons <- readRDS("data_processed/pollen_chronologies.RDS")

# Isolate full and part sequence dataset_ids
dataset_id_full <- chrons[["dataset_ids"]][["full_sequences"]]
dataset_id_part <- chrons[["dataset_ids"]][["part_sequences"]]
all_dataset_ids <- c(dataset_id_full, dataset_id_part)

# Assign full chronology df to object
pollen_chrons <- chrons[["original_chrons"]]

# Filter out part sequences and assign to object
pollen_part_chrons <-
  pollen_chrons %>% 
  dplyr::filter(dataset_id %in% dataset_id_part)

## Need to determine the top and bottom depths for the sequences in pollen_part

# Make output vectors for loop
dataset_id <- vector()
topdepth <- vector()
botdepth <- vector()

# ID oldest age permissible 
max_age <- 
  10086 +     # Approximate uncal date of Holocene onset
  (4 * 3000)  # Maximum distance between min number of date durations

max_age <- 
  max_age / 100

for (i in unique(pollen_part_chrons$dataset_id)) {
  
  agefile <- pollen_part_chrons %>% 
    dplyr::filter(dataset_id == i) %>% 
    as.data.frame()

  agefile <- agefile %>% 
    dplyr::filter(age_uncalibrated_kyr_bp < max_age)  # filter out samples (4 x 3) + 11.7 calky old

  agefile$diff <- c(diff(agefile$age_uncalibrated_kyr_bp), 10)
  
  agegood <- subset(agefile, age_uncalibrated_kyr_bp > 0.1 & diff < 3)

  row_name <- as.numeric(row.names(agegood))
  result_section <- split(row_name,cummax(c(1,diff(row_name))))
  r1 <- as.vector(rapply(result_section,length,how = "unlist"))
  m <- which(r1==max(r1))
  section <- c(result_section[[m]],result_section[[m]][length(result_section[[m]])]+1)
  agefile_new <- agefile[section,]
  
  topdepth[i] <- agefile_new$depth_cm[1]
  botdepth[i] <- agefile_new$depth_cm[nrow(agefile_new) -1]
  dataset_id[i] <- agefile$dataset_id[1]
  
  
}

top_bot_part_seqs <- 
  as.data.frame(
    cbind(
      dataset_id,
      topdepth,
      botdepth
    )
  ) %>% 
  tidyr::drop_na()

# Save
saveRDS(top_bot_part_seqs,
        "data_processed/top_bottom_depths_part_sequences.RDS")












