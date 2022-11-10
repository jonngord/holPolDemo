# This script processes the radiocarbon SPD data - used as the 
# arch-dates predictor variable in the Holocene diversity ARIMA models   

# Script run on the University of York's research server

library(tidyverse)
library(rcarbon)

# Adapted rcarbon::binsense() (https://cran.r-project.org/web/packages/rcarbon/index.html)
# function to return the cal_bp dates and spd values for each bin size
binsense_adapted <-  function (x, y, h, timeRange, calendar = "BP", binning = "CRA", 
                               raw = F, verbose = T, legend = T, ...) {
  if (!calendar %in% c("BP", "BCAD")) {
    stop("Unknown calendar type")
  }
  if (!binning %in% c("CRA", "calibrated")) {
    stop("binning should be either 'CRA' or 'calibrated'")
  }
  years <- timeRange[1]:timeRange[2]
  xlab <- "Years cal BP"
  coln <- numeric(length = length(h))
  xr <- timeRange
  if (calendar == "BCAD") {
    years <- BPtoBCAD(years)
    xlab <- "Years BC/AD"
    xr <- range(years)
    if (all(xr < 0)) {
      xlab <- "Years BC"
    }
    if (all(xr > 0)) {
      xlab <- "Years AD"
    }
  }
  res <- matrix(NA, nrow = length(years), ncol = length(h))
  craAges <- x$metadata$CRA
  if (length(craAges) != length(y)) {
    stop("x and y have different lengths (each calibrated date in x should have a matching location in y)")
  }
  if (verbose) {
    print("Computing SPDs...")
    flush.console()
    pb <- txtProgressBar(min = 1, max = length(h), style = 3)
  }
  for (b in 1:length(h)) {
    if (verbose) {
      setTxtProgressBar(pb, b)
    }
    if (binning == "CRA") {
      bins <- binPrep(sites = y, ages = craAges, h = h[b])
    }
    if (binning == "calibrated") {
      bins <- binPrep(sites = y, ages = x, h = h[b])
    }
    spdtmp <- spd(x, bins = bins, timeRange = timeRange, 
                  spdnormalised = T, verbose = F, ...)
    date <- spdtmp$grid$calBP
    res[, b] <- spdtmp$grid$PrDens
    coln[b] <- paste("h.", h[b], sep = "")
  }
  
  out <- list(
    res = res,           # RETURN MATRIX OF SPDS WITH DIFFERENT BINSIZES   
    date = date
  )
  
  return(out)
  
}



#### Run on local computer


# Set cores
cores <- 4

# Load data

# Radiocarbon data from Bird et al (2020):

# Download from ------>
# DOI: https://doi.org/10.1038/s41597-022-01118-7
# and add to `data_raw`

# Read in
all_radiocarbon <- read_csv("data_raw/p3k14c_scrubbed_fuzzed.csv")

# Filter age range and drop na rows

set.seed(123) # reproducible subset

vars <- c("Long", "Lat", "SiteName", "Age", "Error")
all_radiocarbon <- all_radiocarbon %>% 
  dplyr::filter(Age < 18000) %>% 
  tidyr::drop_na(any_of(vars)) %>% 
  # in full analysis, all data is used and the SPDs are run on the
  # UoY's research servers. But for this demo, just use a subset
  dplyr::sample_n(5000)               
 

# SPD calculation

# Format radiocarbon dataset
sitenames <- stringr::str_to_lower(all_radiocarbon$SiteName) %>% 
  stringr::str_replace_all("[:punct:]", "") 
sitenames <- gsub('\\s+', '', sitenames)

# Add cal curve to df
all_radiocarbon <- 
  all_radiocarbon %>% 
  dplyr::mutate(cc = ifelse(Lat >= 0, 'intcal20', 'shcal20'))

# Calibrate dates
all_radiocarbon_caldates <- 
  calibrate(x = all_radiocarbon$Age,
            errors = all_radiocarbon$Error,
            calCurves = all_radiocarbon$cc,  
            normalised = TRUE,                                      
            ncores = cores)                                         

# Bin dates within 200 years at each site
all_radiocarbon_bins <- binPrep(sites = sitenames, 
                                ages = all_radiocarbon_caldates,
                                h = 200)

binned_cal_dates <- list(
  all_radiocarbon_caldates = all_radiocarbon_caldates,
  all_radiocarbon_bins = all_radiocarbon_bins,
  sitenames = sitenames
)

# Save binned dates
saveRDS(binned_cal_dates,
        "data_processed/global_binned_cal_radiocarbon_dates.RDS")


# Read in data produced above
binned_cal_dates <-
  readRDS("data_processed/global_binned_cal_radiocarbon_dates.RDS")

all_radiocarbon_caldates <- binned_cal_dates[[1]]
all_radiocarbon_bins <- binned_cal_dates[[2]]
sitenames <- binned_cal_dates[[3]]

# Compute SPD
all_radiocarbon_spd_bins <- spd(all_radiocarbon_caldates,
                                bins = all_radiocarbon_bins,
                                timeRange = c(18000, 0))
saveRDS(all_radiocarbon_spd_bins,
        "data_processed/globalspd_binned.RDS")

# Read in data produced above 
all_radiocarbon_spd_bins <- 
  readRDS("data_processed/globalspd_binned.RDS")
grid <- all_radiocarbon_spd_bins[["grid"]]


# Bin radiocarbon date SPD into 200 year bins
radiocarbon_bins <- grid %>%
  # add time bins
  # 0 - 200 = bin 0,
  # 201 - 400 = bin 1... etc
  dplyr::mutate(bin = floor(calBP / 200) * 200) %>%
  dplyr::group_by(bin) %>% #
  dplyr::summarise(mean_dens = mean(PrDens))

# Calculate change value for spd
radiocarbon_diffs = tibble()
for (i in 90:3) {
  
  older <- radiocarbon_bins[i, "mean_dens"]
  younger <- radiocarbon_bins[i-2, "mean_dens"]
  age <- radiocarbon_bins[i-1, "bin"]
  radiocarbonRoC <- abs(older - younger) # absolute value of the difference between 
  # adjacent bins
  radiocarbon_data <- cbind(age, radiocarbonRoC)
  names(radiocarbon_data)[2] <- "radiocarbon_change"
  radiocarbon_diffs <- rbind(radiocarbon_diffs, radiocarbon_data)
  
}

radiocarbon_diffs <- radiocarbon_diffs %>% 
  dplyr::arrange(bin) %>% 
  dplyr::rename(mean_dens_change = radiocarbon_change)

# Join to radiocarbon_diffs (radiocarbon *change* data)
radiocarbon_bins <- dplyr::full_join(radiocarbon_bins, radiocarbon_diffs,
                                     by = "bin")

# Add leads to radiocarbon raw and radiocarbon change data
radiocarbon_bin_leads <- radiocarbon_bins %>% 
  dplyr::arrange(bin) %>% 
  dplyr::mutate(mean_radiocarbon_spd_lead01 = dplyr::lead(mean_dens, n = 1)) %>%
  dplyr::mutate(mean_radiocarbon_spd_lead02 = dplyr::lead(mean_dens, n = 2)) %>%
  dplyr::mutate(mean_radiocarbon_spd_lead03 = dplyr::lead(mean_dens, n = 3)) %>%
  dplyr::mutate(mean_radiocarbon_spd_lead04 = dplyr::lead(mean_dens, n = 4)) %>%
  dplyr::mutate(mean_radiocarbon_spd_lead05 = dplyr::lead(mean_dens, n = 5)) %>%  
  dplyr::mutate(mean_radiocarbon_spd_lead06 = dplyr::lead(mean_dens, n = 6)) %>%
  dplyr::mutate(mean_radiocarbon_spd_lead07 = dplyr::lead(mean_dens, n = 7)) %>%
  dplyr::mutate(mean_radiocarbon_spd_lead08 = dplyr::lead(mean_dens, n = 8)) %>%
  dplyr::mutate(mean_radiocarbon_spd_lead09 = dplyr::lead(mean_dens, n = 9)) %>%
  dplyr::mutate(mean_radiocarbon_spd_lead10 = dplyr::lead(mean_dens, n = 10)) %>%
  
  dplyr::mutate(mean_radiocarbon_spd_change_lead01 = dplyr::lead(mean_dens_change, n = 1)) %>% 
  dplyr::mutate(mean_radiocarbon_spd_change_lead02 = dplyr::lead(mean_dens_change, n = 2)) %>% 
  dplyr::mutate(mean_radiocarbon_spd_change_lead03 = dplyr::lead(mean_dens_change, n = 3)) %>% 
  dplyr::mutate(mean_radiocarbon_spd_change_lead04 = dplyr::lead(mean_dens_change, n = 4)) %>% 
  dplyr::mutate(mean_radiocarbon_spd_change_lead05 = dplyr::lead(mean_dens_change, n = 5)) %>% 
  dplyr::mutate(mean_radiocarbon_spd_change_lead06 = dplyr::lead(mean_dens_change, n = 6)) %>% 
  dplyr::mutate(mean_radiocarbon_spd_change_lead07 = dplyr::lead(mean_dens_change, n = 7)) %>%  
  dplyr::mutate(mean_radiocarbon_spd_change_lead08 = dplyr::lead(mean_dens_change, n = 8)) %>% 
  dplyr::mutate(mean_radiocarbon_spd_change_lead09 = dplyr::lead(mean_dens_change, n = 9)) %>% 
  dplyr::mutate(mean_radiocarbon_spd_change_lead10 = dplyr::lead(mean_dens_change, n = 10)) %>% 
  
  dplyr::relocate(mean_dens_change, .after = mean_radiocarbon_spd_lead10) %>% 
  dplyr::mutate(bin_age_centre = seq(100, 18100, 200)) %>%
  dplyr::relocate(bin_age_centre, everything()) %>% 
  dplyr::rename(mean_radiocarbon_spd = mean_dens,
                mean_radiocarbon_spd_change = mean_dens_change)

saveRDS(radiocarbon_bin_leads, "data_processed/radiocarbon_data.RDS")



# SPD bin size sensitivity analysis



# Read in data produced above
binned_cal_dates <-
  readRDS("data_processed/global_binned_cal_radiocarbon_dates.RDS")

all_radiocarbon_caldates <- binned_cal_dates[[1]]
all_radiocarbon_bins <- binned_cal_dates[[2]]
sitenames <- binned_cal_dates[[3]]

# Sensitivity analysis -- does the number of bins used change trends in SPDs?
bs <- binsense_adapted(x = all_radiocarbon_caldates,
                       y = sitenames,
                       h = seq(0,500,100),
                       timeRange = c(18000, 0))

# Save binsense data
saveRDS(bs, "data_processed/binsense_data.RDS")






























