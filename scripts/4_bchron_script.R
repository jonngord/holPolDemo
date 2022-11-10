# Set w/d
directory <- "data_processed" 

## Load packages 
packages <- c(
  "tidyverse", "Bchron", "janitor",
  "stringr", "ggpubr", "rlist", 
  "IntCal"
)

for (package in packages) {
  library(package, character.only = TRUE)
  
}

# Set locations for outputted models to be written to
folder_save <- file.path(directory, 
                         paste0("age_depths_bchron"))

if ( ! dir.exists(folder_save)) {
  dir.create(folder_save, showWarnings = FALSE)
}

########################################
######        Data read-in        ######
########################################

# Read in dataset_ids
dataset_ids <- 
  readRDS("data_processed/global_pollen.RDS") %>% 
  dplyr::pull(dataset_id) %>% 
  unique()

# Choose 50  datasets
set.seed(123)
samples <- sample(1:length(dataset_ids), 50)


# Age-depth data from Li et al (2022)
# - https://zenodo.org/record/5793936#.YwdHrXbMI2w

# Download tables that are needed from pangaea
metadata <- read.csv2(
  "data_raw/Table-S1_chronological_control_points_metadata.csv",
  stringsAsFactors = FALSE,
  sep = "\t", dec = ".")

parameter <- read.csv2(
  "data_raw/Table-S3_bacon_parameter_settings.csv",
  stringsAsFactors = FALSE,
  sep = "\t", dec = ".")

AgeDepthPollen <- read.csv2(
  "data_raw/Table-S4_original_chronology_metadata_by_pollen_records.csv",
  stringsAsFactors = FALSE, sep = "\t", dec = ".")

# Read in pollen data 
pollen <- readRDS("data_processed/global_pollen.RDS")

# Additional functions
`%notin%` <- Negate(`%in%`)
numextract <- function(string) {
  str_extract(string, "\\-*\\d+\\.*\\d*")
}

# Define specific Values
delete.chrono <- c("Radiocarbon years BP", "No chronology", NA)
c14 <- c("Carbon-14", "14C", "C14", "14c", "c14", "", "Radiocarbon years BP")
plot.c14.correct <- "age.dating"
plot.AWI <- "AWI.age"


# Run models!
for (i in samples) {
  
  # Isolate dataset_id
  ID <- dataset_ids[i]
  
  ## Make chron tables

  # Chronology table compilation code adapted from Li et al (2022)
  # - https://github.com/LongtermEcology/LegacyAge-1.0/blob/main/R-script_bacon_age-depth_modelling_LegacyAge1_0.R
  
  # creating specific subsets for the ID
  parameter.subset <- parameter[which(parameter[, 1] == ID), ]
  AgeDepthPollen.subset <- AgeDepthPollen[which(AgeDepthPollen$Dataset_ID == ID), ]
  metadata.subset <- metadata[which(metadata$Dataset_ID == ID), ]
  if (is.na(parameter.subset$Reservoir)) {
    parameter.subset$Reservoir <- 0
  }
  if (is.na(parameter.subset$Waterline)) {
    parameter.subset$Waterline <- 0
  }
  calibration <- data.frame(
    as.numeric(metadata.subset$Depth..cm.), 
    as.numeric(metadata.subset$Age_Uncalibrated..kyr.BP.) * 1000, 
    (as.numeric(metadata.subset$Dating.Error_Older..kyr.) + as.numeric(metadata.subset$Dating.Error_Younger..kyr.)) * 500, 
    metadata.subset$Dating_Method, 
    as.numeric(metadata.subset$Age_Calibrated..kyr.BP.) * 1000, 
    as.numeric(metadata.subset$Calibrated.dating_Error..kyr.) * 1000, 
    metadata.subset$Thickness..cm.,
    stringsAsFactors = FALSE)
  names(calibration) <- c("depth", "age", "e.older", "age.type", "cal.age", "cal.age.se", "thickness")
  cal.age <- cal.age.se <- cal.sp <- NULL
  
  # Choosing Calibration Curve
  if (AgeDepthPollen.subset$Latitude..DD.[1] >= 0) {
    cc <- "intcal20"
  }
  if (AgeDepthPollen.subset$Latitude..DD.[1] < 0) {
    cc <- "shcal20"
  }
  if (parameter.subset$Marine) {
    cc <- "marine20"
  }
  
  # set chrons/depths
  chrons <- calibration
  depths <- 
    pollen %>% 
    dplyr::filter(dataset_id == ID) %>% 
    dplyr::distinct(depth_m, .keep_all = TRUE) %>% 
    dplyr::pull(depth_m)
  depths <- depths * 100
  dataset_id <- ID
  
  # If thickness entry == NA, replace with 1
  chrons$thickness[is.na(chrons$thickness)] <- 1
  
  # If error entry <=0, replace with 1
  chrons$e.older[chrons$e.older <= 0] <- 1
  
  # Assign calibration curve
  if(cc == "intcal20") { 
    
    chrons <- 
      chrons %>% 
      dplyr::mutate(
        cal_curve = cc,
        cal_curve = replace(cal_curve, age.type != "Carbon-14", "normal"),
      )  %>%  
      dplyr::mutate(
        cal_curve = replace(cal_curve, age.type == "Dates from prior information", "intcal20")
      )
    
  } else if(cc == "shcal20") {
    
    chrons <- 
      chrons %>% 
      dplyr::mutate(
        cal_curve = cc,
        cal_curve = replace(cal_curve, age.type != "Carbon-14", "normal"),
      )  %>%  
      dplyr::mutate(
        cal_curve = replace(cal_curve, age.type == "Dates from prior information", "shcal20")
      )
    
  } else if(cc == "marine20") {
    
    chrons <- 
      chrons %>% 
      dplyr::mutate(
        cal_curve = cc,
        cal_curve = replace(cal_curve, age.type != "Carbon-14", "normal"),
      )  %>%  
      dplyr::mutate(
        cal_curve = replace(cal_curve, age.type == "Dates from prior information", "marine20")
      )
    
  }
  
  # Set max age
  max_age <- 
    10086 +       # Approximate uncal date of Holocene onset
    ( 4 * 3000 )  # Maximum distance between max number of date durations outside 
  # Holocene
  
  # Filter ages by chron
  chrons <- 
    chrons %>% 
    dplyr::filter(age < max_age) %>% 
    # arrange by depth for bchron
    dplyr::arrange(depth) 
  
  
  
  ########################################
  ######            Model           ######
  ########################################
   
  
  
  # Bchron code adapted from Mottl et al (2021) Rratepol example workflow
  # - https://github.com/HOPE-UIB-BIO/R-Ratepol-package/blob/master/vignettes/workflow-example.Rmd
  
  # Set parameters for Bchron
  i_multiplier <- 1 # change to 5 for full analysis, 1 for speed 
  n_iteration_default <- 10000
  n_burn_default <- 2000
  n_thin_default <- 8
  n_iteration <- n_iteration_default * i_multiplier
  n_burn <- n_burn_default * i_multiplier
  n_thin <- n_thin_default * i_multiplier
  
  # run bchron
  example_01_bchron <- tryCatch(
    {
      Bchron::Bchronology(
        ages = chrons$age,
        ageSds = chrons$e.older,
        positions = chrons$depth,
        calCurves = chrons$cal_curve,
        positionThicknesses = chrons$thickness,
        iterations = n_iteration,
        burn = n_burn,
        thin = n_thin,
        allowOutside = TRUE)
    },
    error = function(e) 
    {
      return(e)
    }
  )
  
  # Predict ages
  age_position <- 
    Bchron:::predict.BchronologyRun(
      example_01_bchron,
      newPositions = depths)
  
  # Predict uncertainties
  age_uncertainties <- 
    age_position %>% 
    as.data.frame() %>% 
    dplyr::mutate_all(., as.integer) %>% 
    as.matrix()
  
  # Reformat
  colnames(age_uncertainties) <- depths
  age_uncertainties_t <- t(age_uncertainties)
  colnames(age_uncertainties_t) <- paste0("draw_", seq(1:1000))
  age_uncertainties_t <- as.data.frame(age_uncertainties_t)
  
  # Use median of age uncertainties as default value of levels
  example_01_level_predicted <- 
    depths %>% 
    as.data.frame() %>% 
    dplyr::rename(depths = 1) %>% 
    dplyr::mutate(
      median_age = apply(
        age_uncertainties, 2,
        stats::quantile,
        probs = 0.5),
      dataset_id = ID)

  # Bind
  bchron_age_depth <- 
    dplyr::bind_cols(
      example_01_level_predicted,
      age_uncertainties_t
    )
  
  # Make list of outputs
  output <- list(dataset_id                     = ID,
                 chron_control_table            = chrons,
                 bchron_age_depth_model         = example_01_bchron,
                 bchron_predicted_median_ages   = example_01_level_predicted,
                 bchron_predicted_uncertainties = age_uncertainties,
                 bchron_age_depth               = bchron_age_depth)
  
  # Set filesave
  file_save <- file.path(folder_save, 
                         paste0("age_depth_", ID, '.RDS'))
  
  # Save
  saveRDS( output , file = file_save )
  
  
}




