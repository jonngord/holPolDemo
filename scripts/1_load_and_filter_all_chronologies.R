# This script loads chronological data for pollen records and filters out records based on 
# the below sampling criteria

source("scripts/load_libraries.R")


## Read in chronologies from Li et al (2022)
# (https://doi.pangaea.de/10.1594/PANGAEA.933132)


# Original chronological information (Table-S1, Li et al)
original_chrons <- 
  readr::read_delim(
    file = "data_raw/Table-S1_chronological_control_points_metadata.csv"
  ) %>% 
  janitor::clean_names()

original_chrons$dating_method %>% unique()

# ID oldest age permissible 
max_age <- 
  10086 +     # Approximate uncal date of Holocene onset
  (4 * 3000)  # Multiply maximum distance (3000y) between max number of date durations outside Hol (4)

max_age <- 
  max_age / 1000

# Filter chronological controls by younger than `max_age`
original_chrons <- 
  original_chrons %>% 
  dplyr::filter(age_uncalibrated_kyr_bp <= max_age)


########################################################################  
#####       Filter pollen sequences by age-control sampling        #####  
########################################################################  


##### Following filtering code adapted from Wang et al (2020)
# DOI: https://doi.org/10.1038/s41597-019-0182-7


################## full core ################## 
####  Choose cores based on age controls & samples
####    criteria: 
####      1. age controls > 5;
####      2. maximum interval < 3000 years;

datasetid <- vector()
for (i in unique(original_chrons$dataset_id)) {
  
  agecontrol <- original_chrons %>% 
    dplyr::filter(dataset_id == i)
  
  number <- nrow(agecontrol)
  
  maxinterval <- max(diff(agecontrol$age_uncalibrated_kyr_bp))
  
  if (is.na(maxinterval)==FALSE)
    if (number >= 5 & maxinterval < 3)
      
      datasetid <- c(datasetid, 
                     agecontrol$dataset_id[1])
  
}


################## part core ################## 
####  Choose cores based on age controls & samples criteria: 
####      age controls that are:
####                  conitnuous;
####                  maximum interval < 3000 years;
####                  in the section pollen samples > 5;

datasetid_part <- vector()
for (i in unique(original_chrons$dataset_id)) {
  
  agecontrol <- original_chrons %>% 
    dplyr::filter(dataset_id == i)
  
  ages <- agecontrol$age_uncalibrated_kyr_bp[which(agecontrol$age_uncalibrated_kyr_bp > 1)]
  agediff <- diff(ages)
  locations <- which(agediff < 3)
  result <- rle(diff(locations))
  
  number <- nrow(agecontrol)
  maxinterval <- max(diff(agecontrol$age_uncalibrated_kyr_bp[which(is.na(agecontrol$age_uncalibrated_kyr_bp)==FALSE)]))
  
  if (maxinterval > 3)
    if (any(result$lengths >= 5 & result$values==1)==TRUE)
      if (number > 5) 
        
        datasetid_part <- c(datasetid_part,
                            agecontrol$dataset_id[1])
  
}

# Join dataset_ids of full and part sequences
dataset_id_keep <- c(datasetid, datasetid_part)
no_seqs <- length(dataset_id_keep)
print(no_seqs)


# Save sequence chronologies that meet inclusion criteria
# Make list of outputs 
chrons <- 
  list(
    dataset_ids = list(
      full_sequences = datasetid,
      part_sequences = datasetid_part
    ),
    original_chrons = original_chrons
  )

saveRDS(chrons, "data_processed/pollen_chronologies.RDS")




