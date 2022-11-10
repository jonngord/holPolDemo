library(tidyverse)
library(patchwork)
library(ggpmisc)
library(tidybayes)

# Read in arima preds
arima_pred_data <- 
  readRDS("data_processed/arima_preds/arima_preds.RDS")

# Extract preds and calculate diffs 
{
  div <- lapply(arima_pred_data, function(div_metric) {
    
    # Extract predicted curves for 1000 resamples
    curve_preds <- lapply(div_metric, function(resample) {
      preds <- resample[["arima_preds"]]
      resample_no <- resample[["resample"]]
      resample_no <- paste0("resample_", resample_no)
      preds$resample <- resample_no
      
      return(preds)
      
    }) %>% 
      dplyr::bind_rows()
    
    # Differences between empirical data and climate models for each resample
    diffs_list_climate <- list()
    for (i in unique(curve_preds$resample)) {
      
      resample_i <- curve_preds %>% 
        dplyr::filter(resample == i)
      
      diff <- 
        (
          resample_i$empirical_data - resample_i$av_climate_only_model
        ) 
      
      diff <- as.data.frame(diff)
      diff$resample <- resample_i$resample
      diff$bin_age_centre <- resample_i$bin_age_centre
      
      diffs_list_climate[[i]] <- diff
      
    }
    
    # Differences between empirical data and human and climate model
    diffs_list_human <- list()
    for (i in unique(curve_preds$resample)) {
      
      resample_i <- curve_preds %>% 
        dplyr::filter(resample == i)
      
      diff <- 
        (
          resample_i$empirical_data - resample_i$av_human_and_climate_model
        )  
      
      
      diff <- as.data.frame(diff)
      diff$resample <- resample_i$resample
      diff$bin_age_centre <- resample_i$bin_age_centre
      
      diffs_list_human[[i]] <- diff
      
    }
    
    outputs <- list(
      curve_preds = curve_preds,
      diffs_list_climate = diffs_list_climate,
      diffs_list_human = diffs_list_human
    )
    
    return(outputs)
    
  })
  
  }

# Set font size 
size = 36

# Plot
{ 
  
  
  # Turnover
  
  # Median preds and empirical
  bc_median_preds_data <-
    div[["turnover"]][["curve_preds"]] %>% 
    dplyr::group_by(bin_age_centre
    ) %>% 
    dplyr::summarise(
      median_human_and_climate_model = median(av_human_and_climate_model),
      median_climate_only_model = median(av_climate_only_model),
      median_human_only_model = median(av_human_only_model),
      median_empirical_data = median(empirical_data)
    ) %>% 
    dplyr::select(!median_human_only_model) %>% 
    tidyr::pivot_longer(!bin_age_centre,
                        names_to = "model",
                        values_to = "bray"
    ) 
  
  # Reorder and rename
 bc_median_preds_data$model <-
   factor(bc_median_preds_data$model,
          levels = c("median_empirical_data", "median_human_and_climate_model", "median_climate_only_model"),
          labels = c('Median \nobserved data', 'Median "human +\nclimate" prediction',  'Median "climate\nonly" prediction')
   )
   
  bc_median_preds_plot <- 
    ggplot(
      data = bc_median_preds_data,
      aes(
        x = bin_age_centre, 
        y = bray,
        colour = model)
    ) + geom_line(
      size = 1.5
    ) + 
    #ylim(0.219, 0.33) +
    scale_x_reverse(limits = c(12000, 0),
                    expand = c(0, 0)) +
    scale_colour_manual(
      name = "",
      values = c("#000000", "#D95F0E", "#31A354"
                 )) + 
    theme_bw() +
    #ylab("Modelled turnover") +
    guides(colour = guide_legend(nrow = 4, byrow = TRUE)) +
    theme(legend.position = "right",
          text = element_text(size = (size)),
          panel.border = element_rect(colour = "black", fill = NA, size = 1), 
          axis.text=element_text(size=size-4),
          legend.text=element_text(size=size-4),
          plot.margin = margin(t = 10, b = 10, r = 10, l = 10, unit = "pt"),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          panel.grid.major.x = element_line(colour="darkgrey"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank()) +
    annotate("text", x = 11000, y = 0.318, label = "F", size = size - 20, fontface = "bold")
  
  
  # Empirical plot
  bc_empirical_interval_plot <- 
    div[["turnover"]][["curve_preds"]] %>% 
    dplyr::select(bin_age_centre, empirical_data) %>%
    group_by(bin_age_centre) %>% 
    median_qi(empirical_data, .width = c(.5, .8, .95)) %>%
    ggplot(aes(x = bin_age_centre, y = empirical_data, ymin = .lower, ymax = .upper)) +
    geom_lineribbon() +
    scale_fill_brewer(palette = "Blues") +
    scale_x_reverse(limits = c(12000, 0),
                    expand = c(0, 0)) +
    #ylim(0.219, 0.33) +
    #ylab("Observed turnover") +
    guides(fill = guide_legend(title = "Bootstrap\n intervals")) +
    theme_bw() +
    theme(legend.title.align = 0.5,
          text = element_text(size = size),
          panel.border = element_rect(colour = "black", fill = NA, size = 1), 
          axis.text=element_text(size=size-4),
          plot.margin = margin(t = 25, b = 10, r = 10, l = 10, unit = "pt"),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          panel.grid.major.x = element_line(colour="darkgrey"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.title = element_text(hjust = 0.5,
                                    vjust = 3)) +
    ggtitle("Turnover") +
    annotate("text", x = 11000, y = 0.318, label = "C", size = size - 20, fontface = "bold")
  
  
  # Cimate diffs 
  bc_preds_diffs_plot_climate_interval <- 
    div[["turnover"]][["diffs_list_climate"]] %>% 
    dplyr::bind_rows() %>% 
    dplyr::select(bin_age_centre, diff)  %>%
    group_by(bin_age_centre) %>% 
    median_qi(diff, .width = c(.5, .8, .95)) %>%
    ggplot(aes(x = bin_age_centre, y = diff, ymin = .lower, ymax = .upper)) +
    geom_lineribbon() +
    scale_fill_brewer(palette = "Greens")+
    geom_hline(yintercept = 0, colour = "red", size = 1) +
    theme(legend.position = "none") +
    #ylab("Turnover climate\n deviation") +
    xlab("") +
    scale_x_reverse(limits = c(12000, 0),
                    expand = c(0, 0)) +
    #ylim(-0.067, 0.067) +
    guides(fill = guide_legend(title = "Bootstrap\n intervals")) +
    theme_bw() +
    theme(text = element_text(size = size),
          panel.border = element_rect(colour = "black", fill = NA, size = 1), 
          axis.text=element_text(size=size-4),
          plot.margin = margin(t = 10, b = 10, r = 10, l = 10, unit = "pt"),
          axis.title.y=element_blank(),
          panel.grid.major.x = element_line(colour="darkgrey"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(margin = margin(t = 10))) +
    annotate("text", x = 11000, y = 0.053, label = "L", size = size - 20, fontface = "bold")
  
  
  # Human diffs 
  bc_preds_diffs_plot_human_interval <- 
    div[["turnover"]][["diffs_list_human"]] %>% 
    dplyr::bind_rows() %>% 
    dplyr::select(bin_age_centre, diff)  %>%
    group_by(bin_age_centre) %>% 
    median_qi(diff, .width = c(.5, .8, .95)) %>%
    ggplot(aes(x = bin_age_centre, y = diff, ymin = .lower, ymax = .upper)) +
    geom_lineribbon() +
    scale_fill_brewer(palette = "YlOrBr")+
    geom_hline(yintercept = 0, colour = "red", size = 1) +
    #ylab("Turnover overall\n deviation") +
    scale_x_reverse(limits = c(12000, 0),
                    expand = c(0, 0)) +
    #ylim(-0.067, 0.067) +
    guides(fill = guide_legend(title = "Bootstrap\n intervals")) +
    theme_bw() +
    theme(legend.title.align = 0.5,
          text = element_text(size = size),
          panel.border = element_rect(colour = "black", fill = NA, size = 1), 
          axis.text=element_text(size=size-4),
          plot.margin = margin(t = 10, b = 10, r = 10, l = 10, unit = "pt"),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          panel.grid.major.x = element_line(colour="darkgrey"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank()) +
    annotate("text", x = 11000, y = 0.053, label = "I", size = size - 20, fontface = "bold")
  
  
  # Richness 
  
  
  # Median preds and empirical
  rich_median_preds_plot <-
    div[["richness"]][["curve_preds"]] %>% 
    dplyr::group_by(bin_age_centre
    ) %>% 
    dplyr::summarise(
      median_human_and_climate_model = median(av_human_and_climate_model),
      median_climate_only_model = median(av_climate_only_model),
      median_human_only_model = median(av_human_only_model),
      median_empirical_data = median(empirical_data)
    ) %>% 
    dplyr::select(!median_human_only_model) %>% 
    tidyr::pivot_longer(
      !bin_age_centre,
      names_to = "Model",
      values_to = "richness"
    ) %>% 
    ggplot(
      aes(
        x = bin_age_centre, 
        y = richness,
        colour = Model)
    ) + geom_line(
      size = 1.5
    ) + 
    scale_x_reverse(limits = c(12000, 0),
                    expand = c(0, 0)) +
    scale_colour_manual(
      name = "",
      values = c("#000000", 
                 "#D95F0E",  
                 "#31A354"),
      breaks = c(
        "median_empirical_data",
        'median_human_and_climate_model',
        'median_climate_only_model'
      ),
      labels = c(
        'Median data',
        'Median "human + climate" prediction',
        'Median "climate only" prediction')
    ) +  
    #ylim(12, 16.8) +
    ylab("Modelled") +
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(size = size),
          panel.border = element_rect(colour = "black", fill = NA, size = 1), 
          axis.text=element_text(size=size-4),
          plot.margin = margin(t = 10, b = 10, r = 10, l = 10, unit = "pt"),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major.x = element_line(colour="darkgrey"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank()) +
    annotate("text", x = 11000, y = 16.3, label = "D", size = size - 20, fontface = "bold")
  
  
  # Empirical plot
  rich_empirical_interval_plot <- 
    div[["richness"]][["curve_preds"]] %>% 
    dplyr::select(bin_age_centre, empirical_data) %>%
    group_by(bin_age_centre) %>% 
    median_qi(empirical_data, .width = c(.5, .8, .95)) %>%
    ggplot(aes(x = bin_age_centre, y = empirical_data, ymin = .lower, ymax = .upper)) +
    geom_lineribbon() +
    scale_fill_brewer(palette = "Blues") +
    scale_x_reverse(limits = c(12000, 0),
                    expand = c(0, 0)) +
    #ylim(12, 16.8) +
    ylab("Observed") +
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(size = size),
          panel.border = element_rect(colour = "black", fill = NA, size = 1), 
          axis.text=element_text(size=size-4),
          plot.margin = margin(t = 25, b = 10, r = 10, l = 10, unit = "pt"),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major.x = element_line(colour="darkgrey"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(), 
          plot.title = element_text(hjust = 0.5,
                                    vjust = 3)) +
    ggtitle("Richness") +
    annotate("text", x = 11000, y = 16.3, label = "A", size = size - 20, fontface = "bold")
  
  
  
  # Climate diffs
  rich_preds_diffs_plot_climate_interval <- 
    div[["richness"]][["diffs_list_climate"]] %>% 
    dplyr::bind_rows() %>% 
    dplyr::select(bin_age_centre, diff)  %>%
    group_by(bin_age_centre) %>% 
    median_qi(diff, .width = c(.5, .8, .95)) %>%
    ggplot(aes(x = bin_age_centre, y = diff, ymin = .lower, ymax = .upper)) +
    geom_lineribbon() +
    scale_fill_brewer(palette = "Greens")+
    geom_hline(yintercept = 0, colour = "red", size = 1) +
    ylab("Climate deviation") +
    xlab("") +
    scale_x_reverse(limits = c(12000, 0),
                    expand = c(0, 0)) +
   # ylim(-4, 4) +
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(size = size),
          panel.border = element_rect(colour = "black", fill = NA, size = 1), 
          axis.text=element_text(size=size-4),
          plot.margin = margin(t = 10, b = 10, r = 10, l = 10, unit = "pt"),
          panel.grid.major.x = element_line(colour="darkgrey"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(margin = margin(t = 10))) +
    annotate("text", x = 11000, y = 3.2, label = "J", size = size - 20, fontface = "bold")
  
  # Human diffs
  rich_preds_diffs_plot_human_interval <- 
    div[["richness"]][["diffs_list_human"]] %>% 
    dplyr::bind_rows() %>% 
    dplyr::select(bin_age_centre, diff)  %>%
    group_by(bin_age_centre) %>% 
    median_qi(diff, .width = c(.5, .8, .95)) %>%
    ggplot(aes(x = bin_age_centre, y = diff, ymin = .lower, ymax = .upper)) +
    geom_lineribbon() +
    scale_fill_brewer(palette = "YlOrBr")+
    geom_hline(yintercept = 0, colour = "red", size = 1) +
    ylab("Overall deviation") +
    scale_x_reverse(limits = c(12000, 0),
                    expand = c(0, 0)) +
    #ylim(-4, 4) +
    guides(fill = guide_legend(title = "Bootstrap\n intervals")) +   
    theme_bw() +
    theme(legend.title.align = 0.5,
          legend.position = "none",
          text = element_text(size = size),
          panel.border = element_rect(colour = "black", fill = NA, size = 1), 
          axis.text=element_text(size=size-4),
          plot.margin = margin(t = 10, b = 10, r = 10, l = 10, unit = "pt"),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major.x = element_line(colour="darkgrey"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank()) +
    annotate("text", x = 11000, y = 3.2, label = "G", size = size - 20, fontface = "bold")
  
  
  # Evenness
  
  # Median preds and empirical
  ev_median_preds_plot <-
    div[["evenness"]][["curve_preds"]] %>% 
    dplyr::group_by(bin_age_centre
    ) %>% 
    dplyr::summarise(
      median_human_and_climate_model = median(av_human_and_climate_model),
      median_climate_only_model = median(av_climate_only_model),
      median_human_only_model = median(av_human_only_model),
      median_empirical_data = median(empirical_data)
    ) %>% 
    dplyr::select(!median_human_only_model) %>% 
    tidyr::pivot_longer(
      !bin_age_centre,
      names_to = "Model",
      values_to = "evenness"
    ) %>% 
    ggplot(
      aes(
        x = bin_age_centre, 
        y = evenness,
        colour = Model)
    ) + geom_line(
      size = 1.5
    ) + 
    scale_x_reverse(limits = c(12000, 0),
                    expand = c(0, 0)) +
    scale_colour_manual(
      name = "",
      values = c("#000000", 
                 "#D95F0E",  
                 "#31A354"),
      breaks = c(
        "median_empirical_data",
        'median_human_and_climate_model',
        'median_climate_only_model'
      ),
      labels = c(
        'Median data',
        'Median "human + \nclimate" prediction',
        'Median "climate \nonly" prediction') 
    ) +  
    #ylim(0.60, 0.71) +
    #ylab("Modelled evenness") +
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(size = size),
          panel.border = element_rect(colour = "black", fill = NA, size = 1), 
          axis.text=element_text(size=size-4),
          plot.margin = margin(t = 10, b = 10, r = 10, l = 10, unit = "pt"),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          panel.grid.major.x = element_line(colour="darkgrey"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank()) +
    annotate("text", x = 11000, y = 0.699, label = "E", size = size - 20, fontface = "bold")
  
  # Empirical data
  ev_empirical_interval_plot <- 
    div[["evenness"]][["curve_preds"]] %>% 
    dplyr::select(bin_age_centre, empirical_data) %>%
    group_by(bin_age_centre) %>% 
    median_qi(empirical_data, .width = c(.5, .8, .95)) %>%
    ggplot(aes(x = bin_age_centre, y = empirical_data, ymin = .lower, ymax = .upper)) +
    geom_lineribbon() +
    scale_fill_brewer(palette = "Blues") +
    scale_x_reverse(limits = c(12000, 0),
                    expand = c(0, 0)) +
    #ylab("Observed evenness") +
    #ylim(0.60, 0.71) +
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(size = size),
          panel.border = element_rect(colour = "black", fill = NA, size = 1), 
          axis.text=element_text(size=size-4),
          plot.margin = margin(t = 25, b = 10, r = 10, l = 10, unit = "pt"),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          panel.grid.major.x = element_line(colour="darkgrey"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.title = element_text(hjust = 0.5,
                                    vjust = 3)) +
    ggtitle("Evenness") +
    annotate("text", x = 11000, y = 0.699, label = "B", size = size - 20, fontface = "bold")
  
  
  # Climate diffs
  ev_preds_diffs_plot_climate_interval <- 
    div[["evenness"]][["diffs_list_climate"]] %>% 
    dplyr::bind_rows() %>% 
    #dplyr::filter(bin_age_centre > 300) %>% 
    dplyr::select(bin_age_centre, diff)  %>%
    group_by(bin_age_centre) %>% 
    median_qi(diff, .width = c(.5, .8, .95)) %>%
    ggplot(aes(x = bin_age_centre, y = diff, ymin = .lower, ymax = .upper)) +
    geom_lineribbon() +
    scale_fill_brewer(palette = "Greens") +
    geom_hline(yintercept = 0, colour = "red", size = 1) +
    #ylab("Evenness climate\n deviation") +
    xlab("Age (yBP)") +
    scale_x_reverse(limits = c(12000, 0),
                    expand = c(0, 0)) +
    #ylim(-0.06, 0.06) +
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(size = size),
          panel.border = element_rect(colour = "black", fill = NA, size = 1), 
          axis.text=element_text(size=size-4),
          plot.margin = margin(t = 10, b = 20, r = 10, l = 10, unit = "pt"),
          axis.title.y=element_blank(),
          panel.grid.major.x = element_line(colour="darkgrey"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.title.x = element_text(vjust = -0.75),
          axis.text.x = element_text(margin = margin(t = 10))) +
    annotate("text", x = 11000, y = 0.047, label = "K", size = size - 20, fontface = "bold")
  
  # Human diffs
  ev_preds_diffs_plot_human_interval <- 
    div[["evenness"]][["diffs_list_human"]] %>% 
    dplyr::bind_rows() %>% 
    dplyr::select(bin_age_centre, diff)  %>%
    group_by(bin_age_centre) %>% 
    median_qi(diff, .width = c(.5, .8, .95)) %>%
    ggplot(aes(x = bin_age_centre, y = diff, ymin = .lower, ymax = .upper)) +
    geom_lineribbon() +
    scale_fill_brewer(palette = "YlOrBr")+
    geom_hline(yintercept = 0, colour = "red", size = 1) +
    #ylab("Evenness overall\n deviation") +
    scale_x_reverse(limits = c(12000, 0),
                    expand = c(0, 0)) +
    #ylim(-0.06, 0.06) +
    guides(fill = guide_legend(title = "Bootstrap\n intervals")) +
    theme_bw() +
    theme(legend.title.align = 0.5,
          legend.position = "none",
          text = element_text(size = size),
          panel.border = element_rect(colour = "black", fill = NA, size = 1), 
          axis.text=element_text(size=size-4),
          plot.margin = margin(t = 10, b = 10, r = 10, l = 10, unit = "pt"),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          panel.grid.major.x = element_line(colour="darkgrey"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank()) +
    annotate("text", x = 11000, y = 0.047, label = "H", size = size - 20, fontface = "bold")
  
  
  # Name plots
  
  turnover <- 
    (bc_empirical_interval_plot / bc_median_preds_plot / bc_preds_diffs_plot_human_interval / bc_preds_diffs_plot_climate_interval) 
  
  richness <- 
    (rich_empirical_interval_plot / rich_median_preds_plot / rich_preds_diffs_plot_human_interval / rich_preds_diffs_plot_climate_interval) 
  
  evenness <-
    (ev_empirical_interval_plot / ev_median_preds_plot / ev_preds_diffs_plot_human_interval / ev_preds_diffs_plot_climate_interval) 
  
  div_plot <- 
    wrap_plots(richness, evenness, turnover) 
  
}

div_plot

# Save
ggsave("figures/fig_02.jpeg",
       plot = div_plot,
       width = 65,
       height = 55,
       units = "cm",
       dpi = 320)





