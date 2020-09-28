# ------------------------------------------------------------------------------
# BioVision Endpoint Assay Initial Data Analysis Script
# September 28, 2020
# TS O'Leary
# ------------------------------------------------------------------------------

# Load packages
require(tidyverse)

# Load data
nadh <- read_csv(here::here("ADH/data/biovision/NADH_std_09212020.csv")) 


x <- nadh %>%
  pivot_longer(cols = -NADH, names_to = "mins", values_to = "abs") %>%
  mutate(mins = as.numeric(str_remove_all(mins, "_min"))) %>%
  group_by(NADH) %>%
  summarize(abs = median(abs))

# Plot the standard curve ---
ggplot(data = x, aes(x = NADH, y = abs)) +
  geom_smooth(method='lm', formula = y ~ x, se = FALSE) +
  geom_point() +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(..eq.label.., ..rr.label.., 
                                          sep = "*\", \"*")),
                        parse = TRUE, 
                        rr.digits = 5,
                        label.x = 0.85, 
                        label.y = "top") +
  theme_classic() +
  labs(y = "Absorbance @ 450 nm", x = "NADH (nmol)") 



