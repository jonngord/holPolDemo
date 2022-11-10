# Load libraries 

# Set directory to source directory
directory <- '.' 

## Load packages 
librarypath <- file.path(directory, "Rlibs")
if (!dir.exists(librarypath)) {
  dir.create(librarypath)
}
.libPaths(c(librarypath, .libPaths()))
packages <- c(
  "tidyverse", 
  "rootSolve",
  "car", 
  "forecast", 
  "boot",
  "janitor",
  "rnaturalearth",
  "neotoma",
  "ggthemes",
  "patchwork",
  "ggthemes",
  "stringr",
  "ggpubr",
  "rlist",
  "rbacon",
  "IntCal",
  "foreach",
  "parallel",
  "doParallel"
)
for (package in packages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, lib = librarypath,
                     repos = 'https://cloud.r-project.org',
                     dependencies = TRUE)
  }
  library(package, character.only = TRUE)
}
