# ------------------------------------------------------------------------------
# BioVision Endpoint Assay Initial Data Analysis Script
# September 28, 2020
# TS O'Leary
# ------------------------------------------------------------------------------

# Load packages
require(tidyverse)

# Load data
nadh <- read_csv(here::here("ADH/data/biovision/NADH_std_09212020.csv"))
#nadh <- read_csv(here::here("ADH/data/biovision/NADH_std_09282020.csv"))


abs <- read_csv(here::here("ADH/data/biovision/adh_abs_09212020.csv"))
abs <- read_csv(here::here("ADH/data/biovision/adh_abs_09282020.csv"))

# Tidy NADH data
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


# Tidy ADH abs data
y <- abs %>%
  pivot_longer(cols = contains("min"), names_to = "mins", values_to = "abs") %>%
  mutate(mins = as.numeric(str_remove_all(mins, "_min")),
         EtOH = as.factor(EtOH)) %>%
  group_by(EtOH, sample_vol, n_flies, genotype, mins) %>%
  summarize(abs = median(abs)) %>%
  filter(n_flies == 5) %>%
  filter(mins != 0) %>%
  filter(sample_vol == 5)

y1 <- y %>%
  mutate(nmol_NADH = (abs - 0.123) / 0.0683) %>%
  mutate(deltaAbs = abs - abs[mins == 10],
         deltaNADH = nmol_NADH - nmol_NADH[mins == 10],
         deltaTime = mins - 10) %>%
  filter(deltaTime != 0) %>%
  mutate(mU_mL = (deltaNADH/(deltaTime*(sample_vol/1000)))*9)
  
  
  pivot_wider(names_from = mins, values_from = abs)
  summarize()

# Plot the absorbance over time
ggplot(data = y, 
       aes(x = mins, y = abs, 
                     color = EtOH)) +
  geom_smooth(method='lm', formula = y ~ x, se = FALSE) +
  geom_point() +
  theme_classic() +
  labs(y = "Absorbance @ 450 nm", x = "Time (mins)") 


