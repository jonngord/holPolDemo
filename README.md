# Nine millennia of human-associated plant diversity increases

# Authors
Jonathan D. Gordon*, Brennen Fagan, Nicky Milner, Chris D. Thomas. 

*Corresponding author: jdg548@york.ac.uk

This work was funded by a Leverhulme Trust Research Centre - The Leverhulme Centre for 
Anthropocene Biodiversity. Any queries should be sent to the corresponding author.

# Project description
The majority of the analysis scripts for the manuscript 'Nine millennia of human-associated plant diversity increases' 
were run on the  University of York's High Performance Computing cluster, Viking, or the University of York's research servers. 
Without using computing resources such as these, the computational time required to run these analyses is intractable. 
This project (holPol_demo) contains the same analysis scripts that were used to generate the results from the manuscript, 
but they have been tailored to run on a random subset of sites, with a significantly reduced number of resamples (the full analyses 
involve > 1000 pollen datasets, each resampled 1000 times). These scripts will run on your local machine in loops, rather than sent 
to remote machines in parallel. Depending on the number of datasets included (default dataset *n* = 50, but you can change this as desired), 
running all scripts from start to finish will take a minimum of 5 hours.  

**Please note** - Given the subset of sites that this demo project analyses, the outputted results will likely be very different
to those contained in the manuscript. *This demo project therefore aims to highlight the method, rather than reproduce the results*. 
The greater the number of datasets you include, the more similar the results will look to the manuscript's.
To reproduce the results **exactly** from the manuscript, you will need to rereun the scripts from the analysis **in full** using a HPC cluster.

## Data 

Pollen count data: https://doi.org/10.1594/PANGAEA.929773,
Chronological data for pollen cores: https://doi.org/10.1594/PANGAEA.933132,
Temperature data: https://doi.org/10.1038/s41586-021-03984-4,
Precipitation data:  https://doi.org/10.1038/s41597-020-0552-1,
Arch-dates; global archaeologically attested radiocarbon dates: 
https://doi.org/10.1038/s41597-022-01118-7.,


All analyses are run in R (R: A language and environment for statistical 
computing. R Foundation for Statistical Computing, Vienna, Austria.URL
https://www.R-project.org/.)


Download data files as promted and save them to 'data_raw`. 
Run the scripts in the following order:


## Code

## Setup

0_project_setup.R (this script creates the file structure required to run analyses)

## Pollen processing

1_load_and_filter_all_chronologies.R (this script loads chronological data for pollen records from Li et al (2022) and filters out records based on sampling criteria)   

2_load_all_pollen_regions.R (this script loads all pollen count data from Herzschuh et al (2021) and filters out sites that didn't make the filter from load_and_filter_all_chronologies.R)

3_id_top_bottom-depths_part_seqs.R (this script identifies the top and bottom depths at which to slice records that only partially met the chronological criteria)   

## Age-depth modelling

4_bchron_script.R*** (this script computes the age-depth models for all pollen records that made it through the filters of previous scripts)

## Pollen resampling

5_viking_resample_global_pollen.R*** (this script resamples 150 pollen grains from each sample in each pollen record)

6_join_resamples_and_age_model_draws.R** (this script joins a draw from the posterior distribution of each of the pollen record's age-depth model with one of the 1000 resampled global pollen datasets)

## Raw diversity calculations

7_viking_diversity_analysis_resamples.R*** (this script computes the diversity (richness, evenness, turnover [as measured by Bray-Curtis]) of each of the 1000 resampled pollen datasets)

## Process palaeo data

8_precipitation_data_prep.R (this script processes the precipitation data used as a predictor variable in the Holocene diversity ARIMA models)

9_temperature_data_prep.R (this script processes the temperature data used as a predictor variable in the Holocene diversity ARIMA models)

10_spd_data_prep.R** (this script processes the radiocarbon SPD data used as the arch-dates predictor variable in the Holocene diversity ARIMA models)

## ARIMA modelling

11_viking_richness_global_arima.R*** (this script models Holocene pollen-type richness and makes counterfactual predictions from the chosen model)

12_viking_evenness_global_arima.R*** (this script models Holocene pollen-type evenness and makes counterfactual predictions from the chosen model)

13_viking_bc_global_arima.R*** (this script models Holocene pollen-type richness and makes counterfactual predictions from the chosen model)

14_extract_arima_preds.R (this script reads the files outputted by by the ARIMA modelling scripts and saves them to a list)

## Main text figures

15_fig_1R (this script plots Fig. 1 from the manuscript)
 
16_fig_2.R (this script plots Fig. 2 from the manuscript)

17_fig_3.R (this script plots Fig. 3 from the manuscript)

## Supplementary figures and tables

18_supplementary_figs.R (this script plots the supplementary figures from the manuscript)











