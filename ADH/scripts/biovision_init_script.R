# ------------------------------------------------------------------------------
# BioVision Endpoint Assay Initial Data Analysis Script
# September 28, 2020
# TS O'Leary
# ------------------------------------------------------------------------------

# Load packages
require(tidyverse)

# Source functions
source(here::here("ADH/scripts/adh_batch_functions.R"))

# Load and tidy data
setwd(here::here("ADH/data/biovision"))
nadh <- read_tidy_nadh_csv("NADH_std_10012020.csv")
abs <- read_tidy_adh_abs_csv("adh_abs_10012020.csv")

# Get the means of the NADH data
nadh_med <- nadh %>%
  group_by(NADH) %>%
  summarize(abs = median(abs))

# Plot the standard curve 
plot_std_curve(nadh_med)

# Get the slopes and intercepts of the standard curve
std_curve_lm <- std_curve(nadh_med)

# Get specific ADH data.
abs_0_f <- abs %>%
  mutate(EtOH = as.factor(EtOH)) %>%
  group_by(EtOH, temp, sample_vol, n_flies, genotype, mins) %>%
  summarize(abs = median(abs)) %>%
  filter(mins %in% c(0, 60))

# Plot the absorbance over time
plot_sample_abs_time(abs_0_f)

# Calculate ADH activity
t <- calc_adh_activity(abs_0_f, std_curve_lm, t_0 = 0)

# Plot the MM Curve
plot_mm_curve_biovision(t %>%
                          filter(EtOH != 0))

calc_Km_Vm_biovision(t)
