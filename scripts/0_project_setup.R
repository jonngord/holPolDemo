# This script creates the required file structure

# Step 1:

# - Set working directory to desired location/make project 



# Step 2

# - Make directories 
# Set working directory
fol_data_raw <- file.path("data_raw")
fol_data_processed <- file.path("data_processed")

if ( ! dir.exists(fol_data_raw)) {
  dir.create(fol_data_raw, showWarnings = FALSE)
}

if ( ! dir.exists(fol_data_processed)) {
  dir.create(fol_data_processed, showWarnings = FALSE)
}




# Step 3

# Download raw data and add to folder `data_raw' when prompted

# Pollen count data: https://doi.org/10.1594/PANGAEA.929773,
# Chronological data for pollen cores: https://doi.org/10.1594/PANGAEA.933132,
# Temperature data: https://doi.org/10.1038/s41586-021-03984-4,
# Precipitation data:  https://doi.org/10.1038/s41597-020-0552-1,
# Arch-dates; global archaeologically attested radiocarbon dates: 
#  https://doi.org/10.1038/s41597-022-01118-7.,
