# This script loads all pollen count data from Herzschuh et al, 2021 and filters pollen records 
# depending on the chronological filtering step in load_and_filter_all_chronologies.R and the number 
# of pollen grains in each sample

source("scripts/load_libraries.R")

## Read in pollen data from Herzschuh-etal_2021
# DOI: 

# Africa
data_africa <- 
  readr::read_delim(
    file = "data_raw/pollen/Africa_pollen_datasets/pollen_counts_africa.csv"
)

# Asia
data_asia <- 
  readr::read_delim(
    file = "data_raw/pollen/Asia_pollen_datasets/pollen_counts_asia.txt"
  )

# Europe
data_europe <- 
  readr::read_delim(
    file = "data_raw/pollen/Europe_pollen_datasets/pollen_counts_europe.csv"
  )

# Indo-Pacific
data_indopacific <- 
  readr::read_delim(
    file = "data_raw/pollen/Indopacific_pollen_datasets/pollen_counts_indopacific.csv"
  )

# Latin-America
data_latin_america <- 
  readr::read_delim(
    file = "data_raw/pollen/Latin_America_pollen_datasets/pollen_counts_south_america.csv"
  )

# North America East
data_north_america_east <- 
  readr::read_delim(
    file = "data_raw/pollen/North_America_East_pollen_datasets/pollen_counts_north_america_east.csv"
  )

# North America East
data_north_america_west <- 
  readr::read_delim(
    file = "data_raw/pollen/North_America_west_pollen_datasets/pollen_counts_north_america_west.csv"
  )

# Combine into one all-regions df
all_regions <-
  dplyr::bind_rows(data_africa,
                   data_asia,
                   data_europe,
                   data_indopacific,
                   data_latin_america,
                   data_north_america_east,
                   data_north_america_west) %>% 
  janitor::clean_names()

total_n_sites <- length(unique(all_regions$dataset_id)) # 2676
total_n_samples <- nrow(all_regions)                    # 148431

# Read in pollen chronology .RDS
chrons <- readRDS("data_processed/pollen_chronologies.RDS")

# Isolate full and part sequence dataset_ids that passed the chron filtering stage
dataset_id_full <- chrons[["dataset_ids"]][["full_sequences"]]
dataset_id_part <- chrons[["dataset_ids"]][["part_sequences"]]
all_dataset_ids <- c(dataset_id_full, dataset_id_part)

# Filter sites that made the chronological filter
all_regions_filtered <- 
  all_regions %>% 
  dplyr::filter(dataset_id %in% all_dataset_ids)

n_sites <- length(unique(all_regions_filtered$dataset_id))
n_samples <- nrow(all_regions_filtered)

# Filter out sites with <150 pollen grains per sample
# Metadata 
sequence_meta <- 
  all_regions_filtered %>% 
  dplyr::select(event:depth_m)

# Pollen data
sequence_counts <- 
  all_regions_filtered %>% 
  dplyr::select(-c(event:depth_m)) %>% 
  as.data.frame()

# Change NAs to zeros 
sequence_counts[is.na(sequence_counts)] <- 0

# Caclulate row sums and remove all rows with < 150 grains
sequence_counts$row_sums <- rowSums(sequence_counts)
sequence_counts <- 
  sequence_counts %>% 
  dplyr::bind_cols(sequence_meta) %>% 
  dplyr::filter(row_sums >= 150) %>% 
  dplyr::select(- row_sums)

global_pollen <- 
  sequence_counts %>% 
  relocate(event:depth_m, dplyr::everything())

 saveRDS(object = global_pollen,
         file = "data_processed/global_pollen.RDS")



  