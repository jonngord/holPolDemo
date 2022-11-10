library(tidyverse)
library(dagitty)
library(glmmTMB)
library(boot)
library(hrbrthemes)

# Fig S1 and Fig S2
# - Flowcharts 

source("scripts/supplementary_flowcharts.R")

# Save outputs manually in R `Plots` window

# Fig S3
# - DAG (powerpoint)


# Fig S4 
# - LAGS
source("scripts/diversity_analyses/drivers/quantify_lags.R")

# Save lags plot
ggsave("figures/model_lags.jpeg",
       plot = model_lags,
       height = 20,
       width = 20,
       units = "cm",
       dpi = 320)

# Fig. S5
# - Beta regression against time interval plot

### Read in data and isolate model and all_paired_diffs for 1 of the resamples
resample_1_turnover <- 
  readRDS("data_processed/diversity_results/turnover_results/resample_1_turnover.RDS")
model <- resample_1_turnover[["model"]]
all_paired_diffs <- resample_1_turnover[["all_paired_diffs"]]

### Predict from model the expected bc per time interval
preds_population <- 
  data.frame(
    time_diff = all_paired_diffs$time_diff
  )

preds_pop <- predict(model, 
                     newdata = preds_population)
preds_pop_inv <- plogis(preds_pop) 
preds_population$preds_population_inverse <- preds_pop_inv

### Plot
p <- 
  ggplot(all_paired_diffs, aes(x = time_diff, y = bc_diff.scaled))  +
  # geom_point(alpha = 0.02, colour = "cyan4") +
  geom_bin_2d() +
  xlab("Time Difference") + 
  ylab("Bray-Curtis") +
  theme_bw()

p <- 
  p + geom_line(preds_population, 
                mapping = aes(x = time_diff, 
                              y = preds_population_inverse),
                size = 1.5, 
                colour = "darkorange") +
  scale_x_continuous(trans = 'log10') +
  theme(text = element_text(size = 16))

ggsave("figures/supp_fig_05.jpeg",
       plot = p,
       width = 24,
       height = 12,
       units = "cm",
       dpi = 320)

# Fig S6
# - The effect of varying the bin size on the global Holocene radiocarbon SPD

# Load binsense data
bs <- readRDS("data_processed/binsense_data.RDS")

bs_spd <- bs[["res"]] %>% as.data.frame()
colnames(bs_spd) <- seq(0,500,100)
bs_date <- bs[["date"]] %>% as.data.frame()
binsense_spds <- 
  dplyr::bind_cols(bs_date, bs_spd) %>% 
  dplyr::rename(cal_bp = 1)

binsense_spds_plot <- 
  binsense_spds %>% 
  tidyr::pivot_longer(
    cols = !cal_bp,
    names_to = "bin",
    values_to = "spd"
  ) %>% 
  ggplot(aes(x = cal_bp, y = spd, color = as.integer(bin), group = bin)) +
  geom_line() +
  scale_x_reverse() +
  theme_bw() +
  ylab("Normalised SPD") +
  xlab("Years (calibrated BP)") +
  scale_colour_gradient(name = "bin", 
                        low = "yellow", high = "magenta")

rib <- binsense_spds %>% dplyr::select(cal_bp, `100`, `200`, `500`)

binsense_long <- binsense_spds %>% 
  tidyr::pivot_longer(
    cols = !cal_bp,
    names_to = "Bin Size\n(years)",
    values_to = "spd"
  ) %>% 
  dplyr::filter(`Bin Size\n(years)` %in% c(100, 200, 500))

binsense_plot <-
  ggplot() +
  geom_line(data = binsense_long,
            mapping = aes(x = cal_bp, y = spd, color = `Bin Size\n(years)`), alpha = 0.8) +
  geom_ribbon(data = rib,
              mapping = aes(x = cal_bp, ymin = `100`, ymax = `500`),
              alpha = 0.2) +
  scale_x_reverse(limits = c(11700, 500)) +
  scale_colour_manual(values = c("cyan4", "darkorange", "magenta")) +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  ylab("Normalised SPD") +
  xlab("Years (calibrated BP)")

ggsave("figures/supp_fig_s6.jpeg",
       plot = binsense_plot,
       height = 12,
       width = 24,
       units = "cm",
       dpi = 320)


# Table s2 
# Save median lags table
lags_medians %>% 
  dplyr::bind_rows() %>% 
  gt::gt() %>%
  gt::cols_label(
    `Diversity Metric` = md("Diversity metric"),
    Human = md("Arch-dates"),
    Precipitation = md("Precipitation change"), 
    Temperature = md("Temperature change")) %>%
  gt::cols_align(
    align = "center"
  ) %>% 
  gt::gtsave("figures/supp_table_01.png", expand = 2)

# Table s3
source("scripts/proportion_changes_late.R")



