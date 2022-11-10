# Format, clean and group_split pollen data function
clean_data <- function(pollen_dataframe) {
  
  # Assign data to object
  pollen <- pollen_dataframe  
  
  # Delete last, empty row
  pollen <- pollen[-nrow(pollen), ] 
  
  # Delete last, empty row
  pollen <- pollen[-nrow(pollen), ] 
  
  # Set na values to 0
  pollen[is.na(pollen)] <- 0
  
  #Filter pollen recrds by row sums
  pollen <- pollen %>% 
    dplyr::filter(between(age, 0, 11700))
  
  pollen_row_sums <- as_tibble(as.numeric(rowSums(pollen[, -c(1:2)])))
  pollen_row_sums <- pollen_row_sums %>% 
    rename(row_sums = value) %>% 
    as.data.frame()
  
  # Filter out samples with <150 grains
  pollen <- cbind(pollen_row_sums, pollen)
  pollen <- 
    pollen %>% 
    dplyr::filter(row_sums >= 150) 
  
  # Filter out samples with fewer than 5 taxa present
  pollen$taxa_present <- rowSums(pollen[, -c(1:3)] > 0)
  pollen <-
    pollen %>% 
    dplyr::filter(taxa_present > 4) %>% 
    dplyr::select(-c(taxa_present))
  
  # Split 
  pollen <- pollen %>% 
    as_tibble() %>% 
    group_split(dataset_id)
  
  return(pollen)
  
}

# Check records have more than 5 samples


# Bray-Curtis function
run_bray_curtis <- function(pollen_site) {
  
  # Isolate ages and id
  pollen_site_ages <- pollen_site$age
  pollen_site_dataset_id <- pollen_site$dataset_id
  
  # Remove ages from pollen df
  pollen_site <- pollen_site %>% 
    dplyr::select(-c(age, dataset_id, row_sums)) %>% 
    as.data.frame()
  
  # Add ages as row names to pollen df
  rownames(pollen_site) <- pollen_site_ages
  
  # Compute Bray-Curtis
  raw_bc <- vegdist(pollen_site, method = 'bray')
  raw_bc <- as.data.frame(as.matrix(raw_bc))
  
  # Add ages as column 
  raw_bc <- rownames_to_column(raw_bc, var = "time") 
  
  return(raw_bc)
  
  
}

# Compute distance matrix function
distance_computation <- function(raw_bray_curtis_matrix) {
  
  # Work out the distance between all the rows of the data matrix (ttps://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/dist)
  time_diff_mat <- as.matrix(dist(raw_bray_curtis_matrix[,"time"],
                                  diag = TRUE, upper = TRUE))
  
  # Only need differences off the diagonal. So, match matrix sizes by removing the 'time' column first
  time <- raw_bray_curtis_matrix[,1] 
  bc_diff <- raw_bray_curtis_matrix[,-1] 
  
  # Add a row to the time_diff_matrix and check that the matrix columns and rows have the same dimensions, using '=='.
  delta <- (row(time_diff_mat) + 1) == col(time_diff_mat) 
  
  successive_time_diff <- time_diff_mat[delta] # Time interval matrix
  successive_bc_diff <- bc_diff[delta] # Successive bray-curtis values
  
  # Combine differences
  pollen_site_dataset_id <- pollen_site_dataset_id[-1]
  successive_diffs <- data.frame("time" = raw_bray_curtis_matrix[-1,"time"], # drop the first, because we need successive differences
                                 "bc_diff" = successive_bc_diff,
                                 "time_diff" = successive_time_diff,
                                 "dataset_id" = pollen_site_dataset_id)
  
  return(successive_diffs)
  
}

# Bray-Curtis, successive diff function calc
bray_curtis_dissim <- function(pollen_site) {
  
  pollen_site <- dplyr::distinct(pollen_site, age_draw, .keep_all = TRUE)
  
  # Isolate ages and id
  pollen_site_age_draw <- pollen_site$age_draw
  sample_id <- pollen_site$sample_id
  pollen_site_dataset_id <- pollen_site$dataset_id[1]
  pollen_site_meta <- pollen_site[, c(1:7)]
  region <- pollen_site_meta$region[1]
  
  # Calculate midpoint of time for plotting successive diffs 
  age_draw_midpoint <- rollapply(pollen_site_age_draw, 2, mean)
  
  # Add ages as row names to pollen df
  #rownames(pollen_site) <- pollen_site_age_draw
  
  # Remove ages from pollen df
  pollen_site <- pollen_site %>% 
    dplyr::select(-c(1:7)) 
  
  # pollen_site <- pollen_site[rowSums(pollen_site[])>0,]
  
  # Compute Bray-Curtis
  raw_bc <- vegdist(pollen_site, method = 'bray')
  raw_bc <- as.data.frame(as.matrix(raw_bc))
  
  # add column and row names
  colnames(raw_bc) <- pollen_site_age_draw
  raw_bc <- cbind(pollen_site_age_draw, raw_bc) %>% 
    dplyr::rename(time = 1)
  
  # Work out the distance between all the rows of the data matrix 
  # (https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/dist)
  time_diff_mat <- as.matrix(dist(raw_bc[,"time"],
                                  diag = TRUE, upper = TRUE))
  
  # Only need differences off the diagonal. So, match matrix sizes by removing the 'time' column first
  time <- raw_bc[,1] 
  bc_diff <- raw_bc[,-1] 
  
  # Add a row to the time_diff_matrix and check that the matrix columns and rows have the same dimensions, using '=='.
  delta <- (row(time_diff_mat) + 1) == col(time_diff_mat) 
  
  successive_time_diff <- time_diff_mat[delta] # Time interval matrix
  successive_bc_diff <- bc_diff[delta] # Successive bray-curtis values
  
  # Combine differences
  successive_diffs <- data.frame("time" = raw_bc[-1,"time"], # drop the first, because we need successive differences
                                 "bc_diff" = successive_bc_diff,
                                 "time_diff" = successive_time_diff,
                                 "dataset_id" = pollen_site_dataset_id,
                                 "region" = region, 
                                 "age_draw_midpoint" = age_draw_midpoint,
                                 "sample_id" = sample_id[-1])
  
  # Get all pairwise differences under n years for building the model (n is the largest time interval between samples)
  # The maximum distance between any 2 samples allowed is 3000 years, and the actual max distance is 2862.5 years
  lt <- ((time_diff_mat < 3000) * (time_diff_mat != 0)) == 1 
  all_paired_diffs <- data.frame("time_diff" = time_diff_mat[lt],
                                 "bc_diff" = as.matrix(bc_diff)[lt],
                                 "dataset_id" = pollen_site_dataset_id,
                                 "region" = region)
  
  
  return(list(successive_diffs = successive_diffs,
              all_paired_diffs = all_paired_diffs,
              age_draw = pollen_site_age_draw))
  
}


## Modelling functions
# function to transform fraction with 0 and 1 to <1 and >0
transform01 <- function(x) {
  (x * (length(x) - 1) + 0.5) / (length(x))
}

# function to backtransform transform01
backtransform.est <- function(x, n) {
  y <- (x * n - 0.5) / (n - 1)
  return(y)
}

RMSE <- function(x, y) {
  a <- sqrt(sum((y - x)^2) / length(x))
  return(a)
}





