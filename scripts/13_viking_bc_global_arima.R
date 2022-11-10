# This script models Holocene pollen-type richness and makes 
# counterfactual predictions from the chosen model 

packages <- c(
  "tidyverse", "rootSolve",
  "car", "forecast", "boot",
  "qpcR", "foreach", "doParallel"
)

for (package in packages) {
  library(package, character.only = TRUE)
}

# Assign pipe to object
`%>%` <- magrittr::`%>%`

## Set location for outputted data to be written to 
folder_save <- file.path("data_processed/arima_preds")

# Create folder 
if ( ! dir.exists(folder_save)) {
  dir.create(folder_save, showWarnings = FALSE)
}

# Create folder 
if ( ! dir.exists(paste0(folder_save, "/turnover"))) {
  dir.create(paste0(folder_save, "/turnover"), showWarnings = FALSE)
}

folder_save <- ("data_processed/arima_preds/turnover")

###############################################################################
#                             PREPARE DATA                                    #
###############################################################################

# Make list of all 1000 diversity resamples
diversity_data <- as.list(dir("data_processed/diversity_results/turnover_results", 
                              pattern = ".RDS", 
                              full.names = TRUE))
# Setup parallel backend
cores = detectCores()
clust <- parallel::makeCluster((cores/2) - 1)
registerDoParallel(clust)

# All 10 in parallel (3 cores) take ~ 1.5 hours
temp <- foreach(result = 1:length(diversity_data)) %dopar% {
  
  packages <- c(
    "tidyverse", "rootSolve",
    "car", "forecast", "boot",
    "qpcR", "foreach", "doParallel"
  )
  
  for (package in packages) {
    library(package, character.only = TRUE)
  }
  
  ## Subset single resample for analysis
  div <- readRDS(diversity_data[[result]])
  
  
  # Isolate richness
  bray <- div[["successive_diffs"]]
  
  # Clear memory
  rm(div)
  gc()
  
  # Load in temperature data
  temperature <- readRDS("data_processed/global_binned_climate_data.RDS")
  
  # Load in precipitation data
  precipitation <- readRDS("data_processed/global_prec_leads.RDS")
  
  # Load in radiocarbon data 
  radiocarbon <- readRDS("data_processed/radiocarbon_data.RDS")
  # Bind precipitation data to temp data
  clim_spd <- dplyr::full_join(temperature, 
                               precipitation,
                               by = c("bin_age_centre" = "time_bp")) %>% 
    dplyr::distinct()
  
  # Bind climate dataframe to radiocarbon data 
  clim_spd <- dplyr::full_join(clim_spd, 
                               radiocarbon,
                               by = "bin_age_centre") %>% 
    dplyr::distinct()
  
  # Add constant to b/c values to ensure they are between 0,1
  min_resid <- min(bray$resid)
  min_resid <- abs(min_resid) + 0.0001
  bray <- bray %>% 
    dplyr::mutate(fit_resp_constant = resid + min_resid)
  
  # Bin bray-curtis data into 200 year bins
  bray_bins <- bray %>%
    # add time bins
    # 0 - 200 = bin 0, 
    # 201 - 400 = bin 1... etc
    dplyr::mutate(bin = floor(age_draw_midpoint / 200) * 200) %>% 
    dplyr::group_by(bin) %>% 
    dplyr::summarise(bin_mean = mean(fit_resp_constant, na.rm = TRUE)) %>% 
    dplyr::mutate(bin_age_centre = seq(100, 11700, 200)) %>% 
    dplyr::relocate(bin_age_centre, everything()) %>% 
    dplyr::rename(bray_mean_fit = bin_mean) %>% 
    dplyr::mutate(time_bins = seq(200, 11800, 200)) %>% 
    dplyr::relocate(bin_age_centre, bin, time_bins, everything()) %>% 
    tidyr::unite(time_bin, bin:time_bins, sep = "-")
  
  # Join radiocarbon/climate data with bray-curtis data 
  bray_climate_radiocarbon_all <- 
    dplyr::left_join(bray_bins, clim_spd,
                     by = "bin_age_centre") 
  
  head(bray_climate_radiocarbon_all)
  
  ###############################################################################
  #                       Whole Holocene ARIMA                                  #
  ###############################################################################
  
  
  bray_clim <- bray_climate_radiocarbon_all %>% 
    dplyr::arrange(desc(bin_age_centre)) %>% 
    dplyr::slice(- nrow(bray_climate_radiocarbon_all))
  
  # Logit transformation of outcome variable, diversity
  bray_clim$bray_logit <- car::logit(bray_clim$bray_mean_fit)
  
  # With predictor variables
  whole_hol <- as.matrix(
    cbind(bray_clim[, c(4:19,      # gmst
                        31:41,     # prec
                        43:53)])   # C14
  )
  
  whole_hol <- whole_hol %>%
    tibble::as_tibble() %>%
    dplyr::relocate(contains("radiocarbon"), .after = prec_change_lead10) %>%
    as.matrix()
  
  
  ###### Train model on whole Holocene with humans,
  ###### precipitation and temperature in the model
  
  # Column ranges - make a table of all the possible
  # combinations of the ranges of columns:
  # 16 temp, 11 prec., 11 human.
  # All models must include at least temp, prec. and C14 (non-lagged)
  truth_table <- expand.grid(temperature = 1:16, precipitation = 1:11, human = 1:11)
  
  arima_outputs <- list()
  for (i in 1:nrow(truth_table)) {             # all combinations of ranges
    # for (i in 1:100) {                       # test loop
    
    # print(i)
    truth_row <- truth_table[i, ]
    
    if (truth_row$precipitation != 0 & truth_row$human != 0) {
      cols <- c(0 + 1:truth_row$temperature, 16 + 1:truth_row$precipitation, 27 + 1:truth_row$human)
      
    } else if (truth_row$precipitation != 0) {
      cols <- c(0 + 1:truth_row$temperature, 16 + 1:truth_row$precipitation)
      
    } else {
      cols <- (0 + 1:truth_row$temperature)
      
    }
    
    sub_mat <- whole_hol[, cols]
    
    arima_fit_whole <- auto.arima(bray_clim[, "bray_logit"],
                                  xreg = sub_mat,
                                  seasonal = FALSE)
    
    output <- list(
      truth = truth_row,
      model = arima_fit_whole
    )
    
    arima_outputs[[i]] <- output
    
  }
  
  # Rename list elements
  names(arima_outputs) <- paste0("model", "_", 1:length(arima_outputs))
  
  # Extract AICc
  aic_results <- lapply(arima_outputs, function(mod){
    
    aic_val <- mod[["model"]][["aicc"]]
    
  })
  
  # Extract best models by AICc
  aic_results_bind <- t(dplyr::bind_rows(aic_results)) %>%
    as.data.frame()
  
  aic_results_bind$mod <- rownames(aic_results_bind)
  
  aic_results_bind <- aic_results_bind %>%
    dplyr::rename(AIC = 1) %>%
    dplyr::arrange(AIC)
  
  min_aic <- min(aic_results_bind$AIC)
  
  # ID all models that are within 2 AIC of best model
  delta_aics <-
    aic_results_bind %>%
    dplyr::mutate(delta_aic = AIC - min_aic)
  
  best_mod_aics <-
    delta_aics %>%
    dplyr::filter(delta_aic < 2)
  
  # Pull names to loop through
  top_mods <-
    best_mod_aics %>%
    dplyr::pull(mod)
  
  # Calculate Akaike weights for best models
  weights <-
    akaike.weights(best_mod_aics$AIC)$weights %>%
    as.data.frame() %>%
    dplyr::bind_cols(top_mods) %>%
    dplyr::rename(weight = 1,
                  model_name = 2)
  
  
  
  ###############################################################################
  #                                FORECAST                                     #
  ###############################################################################
  
  
  # Forecast for each of the top_mods in turn
  preds <- list()
  for (mod in (top_mods)) {
    
    print(mod)
    
    # Assign best model to convenience object
    best_output <- arima_outputs[[mod]]
    best_mod <- best_output[["model"]]
    best_truth <- best_output[["truth"]]
    
    # Isolate columns used in model fit for preds
    if (best_truth$precipitation != 0 & best_truth$human != 0) {
      best_cols <- c(0 + 1:best_output$truth$temperature, 16 + 1:best_output$truth$precipitation, 27 + 1:best_output$truth$human)
      
    } else if (best_truth$precipitation != 0) {
      best_cols <- c(0 + 1:best_output$truth$temperature, 16 + 1:best_output$truth$precipitation)
      
    } else {
      best_cols <- (0 + 1:best_output$truth$temperature)
      
    }
    
    sub_mat_forecast <- whole_hol[, best_cols]
    model_variables <- colnames(sub_mat_forecast)
    
    # Forecast for whole Holocene, with humans left IN
    fcast_whole <- forecast::forecast(best_mod, xreg = sub_mat_forecast)
    
    fcast_whole <- fcast_whole %>%
      tibble::as_tibble() %>%
      # Convert from logit
      dplyr::mutate(forecast_inv = inv.logit(`Point Forecast`),
                    quant_95 = inv.logit(`Hi 95`),
                    quant_05 = inv.logit(`Lo 95`)) %>%
      # bind with ages
      dplyr::bind_cols(bray_clim$bin_age_centre) %>%
      dplyr::rename(age = 9)
    
    whole_hol_preds <-
      cbind(
        fcast_whole$age,
        fcast_whole$forecast_inv,
        fcast_whole$quant_95,
        fcast_whole$quant_05) %>%
      as.data.frame() %>%
      dplyr::rename(bin_age_centre = V1,
                    bray_forecast = V2,
                    bray_forecast_95 = V3,
                    bray_forecast_05 = V4)
    
    
    # Forecast for whole Holocene with humans set to 0
    sub_mat_no_hum <- sub_mat_forecast %>%
      dplyr::as_tibble() %>%
      dplyr::mutate_at(vars(contains("radiocarbon")), ~ 0) %>%
      as.matrix()
    
    fcast_whole_no_hum <- forecast::forecast(best_mod, xreg = sub_mat_no_hum)
    
    fcast_whole_no_hum <- fcast_whole_no_hum %>%
      tibble::as_tibble() %>%
      # Convert from logit
      dplyr::mutate(forecast_inv = inv.logit(`Point Forecast`),
                    quant_95 = inv.logit(`Hi 95`),
                    quant_05 = inv.logit(`Lo 95`)) %>%
      # bind with ages
      dplyr::bind_cols(bray_clim$bin_age_centre) %>%
      dplyr::rename(age = 9)
    
    whole_hol_preds_no_hum <-
      cbind(
        fcast_whole_no_hum$age,
        fcast_whole_no_hum$forecast_inv,
        fcast_whole_no_hum$quant_95,
        fcast_whole_no_hum$quant_05) %>%
      as.data.frame() %>%
      dplyr::rename(bin_age_centre = V1,
                    bray_forecast = V2,
                    bray_forecast_95 = V3,
                    bray_forecast_05 = V4)
    
    
    # Forecast for whole Holocene with total climate set to 0
    sub_mat_no_clim <- sub_mat_forecast %>%
      dplyr::as_tibble() %>%
      dplyr::mutate_at(vars(contains("gmst")), ~ 0) %>%
      dplyr::mutate_at(vars(contains("prec")), ~ 0) %>%
      as.matrix()
    
    fcast_whole_no_clim <- forecast::forecast(best_mod, xreg = sub_mat_no_clim)
    
    fcast_whole_no_clim <- fcast_whole_no_clim %>%
      tibble::as_tibble() %>%
      # Convert from logit
      dplyr::mutate(forecast_inv = inv.logit(`Point Forecast`),
                    quant_95 = inv.logit(`Hi 95`),
                    quant_05 = inv.logit(`Lo 95`)) %>%
      # bind with ages
      dplyr::bind_cols(bray_clim$bin_age_centre) %>%
      dplyr::rename(age = 9)
    
    whole_hol_preds_no_clim <-
      cbind(
        fcast_whole_no_clim$age,
        fcast_whole_no_clim$forecast_inv,
        fcast_whole_no_clim$quant_95,
        fcast_whole_no_clim$quant_05) %>%
      as.data.frame() %>%
      dplyr::rename(bin_age_centre = V1,
                    bray_forecast = V2,
                    bray_forecast_95 = V3,
                    bray_forecast_05 = V4)
    
    
    # Forecast for whole Holocene with everything set to 0 apart from temperature
    sub_mat_temp_only <- sub_mat_forecast %>%
      dplyr::as_tibble() %>%
      dplyr::mutate_at(vars(contains("prec")), ~ 0) %>%
      dplyr::mutate_at(vars(contains("radiocarbon")), ~ 0) %>%
      as.matrix()
    
    fcast_whole_temp_only <- forecast::forecast(best_mod, xreg = sub_mat_temp_only)
    
    fcast_whole_temp_only <- fcast_whole_temp_only %>%
      tibble::as_tibble() %>%
      # Convert from logit
      dplyr::mutate(forecast_inv = inv.logit(`Point Forecast`),
                    quant_95 = inv.logit(`Hi 95`),
                    quant_05 = inv.logit(`Lo 95`)) %>%
      # bind with ages
      dplyr::bind_cols(bray_clim$bin_age_centre) %>%
      dplyr::rename(age = 9)
    
    whole_hol_preds_temp_only <-
      cbind(
        fcast_whole_temp_only$age,
        fcast_whole_temp_only$forecast_inv,
        fcast_whole_temp_only$quant_95,
        fcast_whole_temp_only$quant_05) %>%
      as.data.frame() %>%
      dplyr::rename(bin_age_centre = V1,
                    bray_forecast = V2,
                    bray_forecast_95 = V3,
                    bray_forecast_05 = V4)
    
    
    # Forecast for whole Holocene with everything set to 0 apart from precipitation
    sub_mat_prec_only <- sub_mat_forecast %>%
      dplyr::as_tibble() %>%
      dplyr::mutate_at(vars(contains("gmst")), ~ 0) %>%
      dplyr::mutate_at(vars(contains("radiocarbon")), ~ 0) %>%
      as.matrix()
    
    fcast_whole_prec_only <- forecast::forecast(best_mod, xreg = sub_mat_prec_only)
    
    fcast_whole_prec_only <- fcast_whole_prec_only %>%
      tibble::as_tibble() %>%
      # Convert from logit
      dplyr::mutate(forecast_inv = inv.logit(`Point Forecast`),
                    quant_95 = inv.logit(`Hi 95`),
                    quant_05 = inv.logit(`Lo 95`)) %>%
      # bind with ages
      dplyr::bind_cols(bray_clim$bin_age_centre) %>%
      dplyr::rename(age = 9)
    
    whole_hol_preds_prec_only <-
      cbind(
        fcast_whole_prec_only$age,
        fcast_whole_prec_only$forecast_inv,
        fcast_whole_prec_only$quant_95,
        fcast_whole_prec_only$quant_05) %>%
      as.data.frame() %>%
      dplyr::rename(bin_age_centre = V1,
                    bray_forecast = V2,
                    bray_forecast_95 = V3,
                    bray_forecast_05 = V4)
    
    # Format empirical data
    empirical <- bray_bins %>%
      dplyr::select(bin_age_centre, bray_mean_fit) %>%
      dplyr::arrange(desc(bin_age_centre)) %>%
      dplyr::filter(bin_age_centre > 100)
    
    
    # Cbind empirical and model preds  together
    intersection_df <- dplyr::bind_cols(
      empirical,
      whole_hol_preds$bray_forecast,
      whole_hol_preds_no_clim$bray_forecast,
      whole_hol_preds_no_hum$bray_forecast,
      whole_hol_preds_prec_only$bray_forecast,
      whole_hol_preds_temp_only$bray_forecast
    ) %>%
      dplyr::rename(
        empirical_data = 2,
        human_and_climate_model = 3,
        human_only_model = 4,
        climate_only_model = 5,
        precipitation_only_model = 6,
        temperature_only_model = 7
      ) %>%
      dplyr::arrange(bin_age_centre) %>%
      dplyr::mutate(model_name = mod)
    
    preds[[mod]] <-  intersection_df
    
  }
  
  
  # Weighted average of model predictions
  preds <-
    preds %>%
    dplyr::bind_rows() %>%
    dplyr::left_join(weights,
                     by = "model_name")
  
  weighted_av_preds <-
    preds %>%
    dplyr::group_by(bin_age_centre) %>%
    summarise(
      av_human_and_climate_model =  sum( human_and_climate_model * weight ),
      av_human_only_model = sum( human_only_model * weight ),
      av_climate_only_model = sum( climate_only_model * weight ),
      av_precipitation_only_model = sum( precipitation_only_model * weight ),
      av_temperature_only_model = sum( temperature_only_model * weight )
    ) %>%
    dplyr::bind_cols(intersection_df$empirical_data) %>%
    dplyr::rename(empirical_data = 7) %>%
    dplyr::relocate(empirical_data, .after = bin_age_centre)
  
  
  
  
  ######### Plot preds ########
  
  
  # Colourblind palette with black:
  cbbPalette <- c("#009E73", "#D55E00","#F0E442", "#CC79A7", "#E69F00","#000000" )
  holocene_diversity_predictions_plot <-
    weighted_av_preds %>%
    pivot_longer(!bin_age_centre,
                 names_to = "model",
                 values_to = "bray") %>%
    ggplot(aes(x = bin_age_centre,
               y = bray,
               colour = model)) +
    geom_line(size = 1.2) +
    scale_x_reverse() +
    theme_bw() +
    xlab("Age (cal. years BP)") +
    ylab("Compositional turnover") +
    scale_colour_manual(values = cbbPalette)
  
  
  # PLEASE NOTE #
  
  # holocene_diversity_predictions_plot (plot just made, above) will begin to look more 
  # the publication plots (Fig. 2) as more datasets_ids are included. This is just a
  # test run with a limited number of i) dataset_ids and ii) resamples, included, so the 
  # results look spikey and quite dissimilar to the publication. 
  
  
  ######### Save outputs #########
  
  
  whole_hol_outputs <-
    list(
      resample = result,
      bray_clim = bray_clim,
      predictor_variables =
        list(
          truth_row = best_truth,
          predictor_variable_names = model_variables
        ),
      aic_results = aic_results_bind,
      arima_model = best_mod,
      holocene_diversity_predictions_plot = holocene_diversity_predictions_plot,
      weighted_av_preds = weighted_av_preds,
      weights = weights
    )
  
  
  
  ######### Save outputs to list #########
  
  both_outputs <-
    list(
      whole_hol_outputs = whole_hol_outputs
    )
  
  
  ## Set filesave
  file_save <- file.path(folder_save, paste0("bray_arima_preds_resample_",
                                             result,
                                             ".RDS"))
  
  ## Save file
  saveRDS(both_outputs,
          file = file_save)
  
  return(result)
  
}







