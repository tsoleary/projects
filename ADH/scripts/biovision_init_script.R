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
x <- nadh %>%
  group_by(NADH) %>%
  summarize(abs = median(abs))

# Plot the standard curve 
plot_std_curve(x)

# Get the slopes and intercepts of the standard curve
std_curve_lm <- std_curve(x)

# Get specific ADH data.
y <- abs %>%
  mutate(EtOH = as.factor(EtOH)) %>%
  group_by(EtOH, temp, sample_vol, n_flies, genotype, mins) %>%
  summarize(abs = median(abs)) %>%
  filter(mins %in% c(0, 60))

# Plot the absorbance over time
plot_sample_abs_time(y)

# Calculate ADH activity
t <- calc_adh_activity(y, std_curve_lm, t_0 = 0)

# Plot the MM Curve
plot_mm_curve(t)
