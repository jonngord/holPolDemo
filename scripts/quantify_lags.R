library(tidyverse)
library(patchwork)

# Read in arima preds
arima_preds <- readRDS("data_processed/arima_preds/arima_preds.RDS")

lags_all <- lapply(arima_preds, function(div_metric){
  
  # For each diversity metric (richness, evenness, turnover)...
  diversity_metric_lags <- lapply(div_metric, function(resample) {
    
    # Extract how many lags were used in the best model for each resample
    lags <- resample[["predictor_variables"]][["truth_row"]]
    lags$resample <- resample[["resample"]]
    
    return(lags)
    
  }) %>% 
    
    # Bind all resamples together
    dplyr::bind_rows() %>% 
    
    # Subtract 1 from each of the datasets because `1` is the raw, 
    # unlagged variable
    dplyr::mutate(Temperature = temperature - 1,
                  Precipitation = precipitation - 1,
                  Human = human - 1) %>% 
    dplyr::select(Temperature:Human) 
  
  # Return lags used in models for each diversity metric
  return(diversity_metric_lags)
  
})

vars <- 
  c(
    `Human` = "Arch-dates",
    `Precipitation` = "Precipitation change",
    `Temperature` = "Temperature change"
  )

# Plot
richness_lags <- 
  lags_all[["richness"]] %>% 
  dplyr::mutate(resample = seq(1, nrow(lags_all[["richness"]]), 1)) %>% 
  tidyr::pivot_longer(!resample, names_to = "variable", values_to = "lags") %>% 
  ggplot(aes(x = lags)) +
  geom_histogram(bins = 40) +
  facet_wrap(~ variable, labeller = as_labeller(vars)) +
  xlab("Number of lags") +
  ylab("Count") +
  ggtitle("Richness Lags") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  xlim(c(-1, 15)) +
  ylim(0, 875)

evenness_lags <-
  lags_all[["evenness"]] %>% 
  dplyr::mutate(resample = seq(1, nrow(lags_all[["evenness"]]), 1)) %>% 
  tidyr::pivot_longer(!resample, names_to = "variable", values_to = "lags") %>% 
  ggplot(aes(x = lags)) +
  geom_histogram(bins = 40) +
  facet_wrap(~ variable, labeller = as_labeller(vars)) +
  xlab("Number of lags") +
  ylab("Count") +
  ggtitle("Evenness Lags") +
  theme_bw()  +  
  theme(text = element_text(size = 16)) +
  xlim(c(-1, 15)) +
  ylim(0, 875)

turnover_lags <- 
  lags_all[["turnover"]] %>% 
  dplyr::mutate(resample = seq(1, nrow(lags_all[["turnover"]]), 1)) %>% 
  tidyr::pivot_longer(!resample, names_to = "variable", values_to = "lags") %>% 
  ggplot(aes(x = lags)) +
  geom_histogram(bins = 40) +
  facet_wrap(~ variable, labeller = as_labeller(vars)) +
  xlab("Number of lags") +
  ylab("Count") +
  ggtitle("Turnover Lags") +
  theme_bw()  +
  theme(text = element_text(size = 16)) +
  xlim(c(-1, 15)) +
  ylim(0, 875)

# Save plots
model_lags <- 
  richness_lags /
  evenness_lags /
  turnover_lags

ggsave("figures/model_lags.jpeg",
       plot = model_lags,
       height = 20,
       width = 24,
       units = "cm",
       dpi = 320)
  
# Medians

# Function to change first letter to caps
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

lags_medians <- list()
for (i in names(lags_all)) {
  
  lags_medians[[i]] <- 
    lags_all[[i]]%>% 
    dplyr::mutate(resample = seq(1, nrow(lags_all[[i]]), 1)) %>% 
    tidyr::pivot_longer(!resample, names_to = "variable", values_to = "lags") %>% 
    group_by(variable) %>% 
    summarise(median_lags = median(lags)) %>% 
    dplyr::mutate(`Diversity Metric` = firstup(i)) %>% 
    tidyr::pivot_wider(names_from = "variable", values_from = "median_lags")
  
}

lags_medians %>% 
  dplyr::bind_rows() %>% 
  gt::gt() %>% 
  gt::gtsave("figures/supp_table_01.png", expand = 10)





















