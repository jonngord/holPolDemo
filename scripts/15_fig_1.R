library(tidyverse)
library(rnaturalearth)
library(ggrepel)
library(ggthemes)
library(patchwork)
library(parallel)
library(doParallel)
library(grid)


####### MAP ####### 

dir.create("figures")

# Make directory of 1000 pollen resamples 
pollen_files <- 
  as.list(
    dir("data_processed/pollen_resample_with_age_draw/",
        pattern = ".RDS",
        full.names = TRUE))


# Setup parallel backend
cores <- parallel::detectCores()
clust <- parallel::makeCluster(cores - 1)
registerDoParallel(clust)

pollen_datasets <- 
  
  foreach(i =  1:length(pollen_files)) %dopar%  {
    
    library(tidyverse)
    
    pollen <- readRDS(pollen_files[[i]])
    
    dataset_ids <- unique(pollen$dataset_id) # dataset_ids in each resample
    
    return(dataset_ids)
    
  }

stopCluster(clust)

# Datset_ids in ALL resamples
pollen_datasets_unlist <- unlist(pollen_datasets)  

# Datasets in study
# - Passed minimum inclusion critera, age-depth modelling, etc 
pollen_datasets_unlist <- unique(pollen_datasets_unlist)  

# Number of datasets in our study across all resamples (some variation between
# resamples given the resampling procedure and subsequent filtering steps)
n_datasets <- length(pollen_datasets_unlist)                

# Load in whole df with attached metadata
all_datasets <- 
  readRDS("data_processed/global_pollen.RDS") 

# Filter by datasets included in the study
pollen_datasets_keep <- 
  all_datasets %>% 
  dplyr::filter(dataset_id %in% pollen_datasets_unlist)

ids <- 
  pollen_datasets_keep %>% 
  dplyr::distinct(dataset_id, .keep_all = T) 




#######  STUDY SITES LIST FOR SUPP MATS #######

# readr::write_csv(ids %>% dplyr::select(1:8), "data_processed/pollen/study_datasets.csv")




####### POLLEN HISTOGRAM #######

# Average numnber of samples per 200 year bin across all resamples

# Setup parallel backend
clust <- parallel::makeCluster(2)
registerDoParallel(clust)

pollen <-
  
  foreach(i =  1:length(pollen_files)) %dopar%  {
    
    library(tidyverse)
    
    pollen <- readRDS(pollen_files[[i]])
    
    resample <- paste0("resample_", i)
    
    sample_n <- pollen %>%
      dplyr::mutate(bin = floor(age_draw / 200) * 200) %>%
      dplyr::group_by(bin) %>%
      dplyr::summarise(sample_n = n()) %>%
      mutate(resample = resample)
    
    return(sample_n)
    
  }

stopCluster(clust)

sample_n <-
  pollen %>%
  dplyr::bind_rows() %>%
  dplyr::group_by(bin) %>%
  summarise(median_sample_n = median(sample_n),
            sd_sample_n = sd(sample_n)) %>%
  dplyr::mutate(bin_age_centre = seq(100, 11700, 200))




####### PALAEO DATA #######



# # Read in palaeo data
global_binned_climate_data <- readRDS("data_processed/global_binned_climate_data.RDS")
global_prec_leads <- readRDS("data_processed/global_prec_leads.RDS")
radiocarbon_data <- readRDS("data_processed/radiocarbon_data.RDS")

# Make df containing all variable for plotting
palaeo_data <-
  tibble::tibble(
    age = global_binned_climate_data$bin_age_centre[1:58],
    temp_change = global_binned_climate_data$gmst_change[1:58],
    prec_change = global_prec_leads$prec_change[1:58],
    radiocarbon = radiocarbon_data$mean_radiocarbon_spd[2:59]
  ) %>%
  dplyr::filter(age >= 300)



####### PALAEO DATA  PLOTS #######



# Set font sizes
size <- 14
size_element <- 18


# Plots

# Pollen histogram
pollen_hist_global <-
  ggplot(sample_n %>% 
           dplyr::filter(bin_age_centre >= 300), 
         aes(x = bin_age_centre,
             y = median_sample_n)) + 
  geom_bar(stat = "identity",
           colour = "darkgrey",
           alpha = 0.9) +
  ylab("Median pollen\nsample count") +
  xlab("Age (yBP)") +
  theme_bw() +
  scale_x_reverse(limits = c(12000, -2500),
                  expand = c(0, 0),
                  breaks = c(12000, 9000, 6000, 3000, 0)) +
  scale_y_continuous(limits = c(0, 1700),
                     expand = c(0, 100)) +
  coord_flip() +
  theme(
    text = element_text(size = size),
    plot.margin = margin(t = 0.5, b = -50, r = 0, l = 0, unit = "cm"),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x=element_text(size = size_element), 
    axis.text.y=element_text(size = size_element),
    panel.grid.major.y = element_line(colour="darkgrey"),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.ticks = element_line(colour = "darkgrey"),
    axis.title.x = element_text(size = size_element),
    axis.title.y = element_text(size = size_element))

pollen_hist_global <- 
  pollen_hist_global + annotate("text", x = -1000, y = 0, label = "B", size = size-7, fontface = "bold")

# Temp change 
temp_plot <-
  ggplot(
    data = palaeo_data,
    mapping = aes(
      x = age,
      y = temp_change)
  ) + 
  geom_line(size = 1, colour = "black") +
  theme_bw() +
  coord_flip() +
  ylab("Absolute\n temperature \nchange (°C)") +
  xlab("") +
  scale_x_reverse(limits = c(12000, -2500),
                  expand = c(0, 0),
                  breaks = c(12000, 9000, 6000, 3000, 0)) +
  scale_y_continuous(limits = c(-0.02, 0.85)) +
  coord_flip() +
  theme(
    text = element_text(size = size),
    plot.margin = margin(t = 0.5, b = 0, r = 0, l = 0, unit = "cm"),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.text.x=element_text(size = size_element), 
    panel.grid.major.y = element_line(colour="darkgrey"),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.ticks = element_line(colour = "darkgrey"),
    axis.line.y = element_line(color="darkgrey", size = 1, linetype = 2),
    axis.title.x = element_text(size = size_element))

temp_plot <-
    temp_plot + annotate("text", x = - 1000, y = 0, label = "C", size = size-7, fontface = "bold")

# Prec change
prec_plot <- 
  ggplot(
    data = palaeo_data,
    mapping = aes(
      x = age,
      y = prec_change)
  ) + 
  geom_line(size = 1, colour = "black") +
  scale_x_reverse(limits = c(12000, -2500),
                  expand = c(0, 0),
                  breaks = c(12000, 9000, 6000, 3000, 0)) +
  theme_bw() +
  coord_flip() +
  scale_y_continuous(limits = c(-0.8, 11), 
                     expand = c(0, 0.2),
                     breaks = seq(0, 10, 2)) +
  ylab("Absolute\n precipitation\n change (mm)") +
  xlab("") +
  theme(
    text = element_text(size = size),
    plot.margin = margin(t = 0.5, b = 0, r = 0, l = 0, unit = "cm"),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.text.x=element_text(size = size_element), 
    panel.grid.major.y = element_line(colour="darkgrey"),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.ticks = element_line(colour = "darkgrey"),
    axis.line.y = element_line(color="darkgrey", size = 1, linetype = 2),
    axis.title.x = element_text(size = size_element)) 


prec_plot <-
  prec_plot + annotate("text", x = - 1000, y = 0, label = "D", size = size-7, fontface = "bold")

# Arch-dates
radiocarbon_plot <- 
  ggplot(
    data = palaeo_data,
    mapping = aes(
      x = age,
      y = radiocarbon)
  ) + 
  geom_bar(stat = "identity", 
           width = 200,
           colour = "darkgrey",
           alpha = 0.9) +
  scale_x_reverse(limits = c(12000, -2500),
                  expand = c(0, 0),
                  breaks = c(12000, 9000, 6000, 3000, 0)) +
  scale_y_continuous(limits = c(-1, 18),
                     breaks = seq(0, 16, 4),
                     expand = c(0, 0.2)) +
  coord_flip() +
  theme_bw() +
  ylab("Arch-dates") +
  xlab("") +
  theme(
    text = element_text(size = size),
    plot.margin = margin(t = 0.5, b = 0, r = 0, l = 0, unit = "cm"),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.text.x=element_text(size = size_element), 
    panel.grid.major.y = element_line(colour="darkgrey"),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.ticks = element_line(colour = "darkgrey"),
    axis.line.y = element_line(color="darkgrey", size = 1, linetype = 2),
    axis.title.x = element_text(size = size_element))

radiocarbon_plot <-
  radiocarbon_plot + annotate("text", x = - 1000, y = 0, label = "E", size = size- 7, fontface = "bold")

# Map
world <- ne_countries(scale = "medium", returnclass = "sf")

map <- 
  ggplot(data = world) +
  geom_sf(colour = NA,
          fill = "darkgray",
          alpha = 0.5) +
  geom_point(data = ids, 
             aes(x = longitude_dd, 
                 y = latitude_dd),
             colour = "#FF7518",
             alpha = 1, 
             size = 2.25) + 
  #colour = region)) +
  coord_sf(ylim = c(-63, 90),
           xlim = c(-185, 185),
           expand = F) +
  xlab("") +
  ylab("") +
  theme_map() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    plot.margin = margin(t = 0, b = 0, r = 0, l = 0, unit = "cm"),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )

map <-
  map + annotate("text", x = -175, y = 82, label = "A", size = size - 7, fontface = "bold")




#  Composite plot with patchwork


# Make labels
df <- data.frame()
labs <- ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 100) + theme_void() +
  annotate("text", x = 0, y = 85, label = "1850 CE", hjust = 0, size = 6) +
  annotate("text", x = 3, y = 0, label = "P/H", hjust = 0, size = 6) 

palaeo <- (pollen_hist_global | temp_plot | prec_plot | radiocarbon_plot | labs) +
  plot_layout(nrow = 1,
              ncol = 5, 
              width = c(3, 3, 3, 3, 1.3)) +
  theme(plot.margin = margin(t = 0, b = 0, r = 0, l = 0, unit = "cm"))

both <- 
  map / palaeo  +
  plot_layout(nrow = 2,
              ncol = 1,
              height = c(7, 3), 
              width = c(9, 9))

# Save
ggsave("figures/fig_01.jpeg",
       plot = both, 
       units = "cm")





