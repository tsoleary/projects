# ------------------------------------------------------------------------------
# ADH assay batch processing functions
# June 10, 2020
# TS O'Leary
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Function: read_tidy_spec_csv
# Description: uses read_csv and tidys the data
# Inputs: .csv file 
# Outputs: data.frame

require(tidyverse)

read_tidy_spec_csv <- function(file) {
  
  read_csv(file) %>%
    pivot_longer(-c("cycle", "time", "temp"), 
                 names_to = "conc", 
                 values_to = "abs") %>%
    mutate(conc = as.numeric(str_replace_all(conc, "_.", ""))) %>%
    mutate(file = file,
           temp_group = case_when(temp > 17.5 & temp < 19.5 ~ "18.5",
                                  temp > 22.5 & temp < 24.5 ~ "23.5",
                                  temp > 27.5 & temp < 29.5 ~ "28.5",
                                  temp > 32.5 & temp < 34.5 ~ "33.5",
                                  temp > 37.5 & temp < 39.5 ~ "38.5"))
  
} 
# End function -----------------------------------------------------------------


# ------------------------------------------------------------------------------
# Function: read_bind_tidy_csv
# Description: Read all the spec csvs in a directory
# Inputs: dir is the directory with csv files to load and bind
# Outputs: large binded tibble

require(tidyverse)

read_bind_tidy_csv <- function(dir) {
  setwd(dir)
  files <- list.files()
  dfs <- vector(mode = "list", length = length(files))
  
  for (i in seq_along(files)){
    dfs[[i]] <- read_tidy_spec_csv(files[i])
  }
  
  bind_rows(dfs)
} 
# End function -----------------------------------------------------------------


# ------------------------------------------------------------------------------
# Function: calc_velocity
# Description: Linear regression to determine velocity of the reaction in 
#              absorbance over time at each conc for each temp_group
# Inputs: tibble from read_bind_tidy_csv() with conc, temp_group, 
#         abs, & time columns
# Outputs: tibble with velocity estimates

require(tidyverse)
require(broom)

calc_velocity <- function(data) {
  data %>%
    group_by(conc, temp_group) %>%
    nest() %>%
    mutate(
      fit = map(data, ~ lm(abs ~ time, data = .x)),
      tidied = map(fit, tidy)
    ) %>% 
    unnest(tidied) %>%
    filter(term == "time" & conc != "0" & estimate > 0) %>%
    select(temp_group, conc, estimate) %>%
    rename(velocity = estimate)
} 
# End function -----------------------------------------------------------------


# ------------------------------------------------------------------------------
# Function: calc_Km_Vm
# Description: Calculate Vmax and Km for Michaelis–Menten curve
# Inputs: tibble with velocity estimates for each conc and temp_group
# Outputs: tibble with estimates for Vmax and Km

require(tidyverse)
require(broom)

calc_Km_Vm <- function(data) {
  data %>%
    filter(conc != 0 & velocity > 0) %>%
    group_by(temp_group) %>%
    nest() %>%
    mutate(
      fit = map(data, ~ nls(formula = 
                              velocity ~ SSmicmen(conc, Vmax, Km), data = .), 
                data = .x),
      tidied = map(fit, tidy)
    ) %>%
    unnest(tidied)
} 
# End function -----------------------------------------------------------------


# ------------------------------------------------------------------------------
# Function: plot_velo
# Description: Plot the absorbance over time for each temperature group at each
#              concentration
# Inputs: tibble from read_bind_tidy_csv function
# Outputs: facet_wrap'd plot

plot_velo <- function(data) {
  ggplot(data, aes(x = time, 
                   y = abs, 
                   color = as.factor(conc))) +
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap(~ temp_group) +
    labs(x = "Time (s)", y = "Absorbance @ 340 nm") +
    scale_color_discrete(name = "Ethanol (mM)") +
    geom_point() +
    theme_classic(base_size = 16)
}
# End function -----------------------------------------------------------------


# ------------------------------------------------------------------------------
# Function: plot_mm_curve
# Description: Plot the Michaelis-Menten curves
# Inputs: tibble with conc and velocity columns
# Outputs: plot

require(tidyverse)

plot_mm_curve <- function(data, 
                          palette = c("skyblue", "blue", "purple", "red", "darkred")) {
  
  ggplot(data, aes(x = conc, y = velocity, fill = temp_group)) +
    geom_point(pch = 21, size = 3, color = "grey50", alpha = 0.75) +
    geom_smooth(method = "nls", 
                formula = y ~ SSmicmen(x, Vmax, Km),
                data = data,
                se = FALSE,
                aes(color = temp_group), 
                show.legend = FALSE) +
    scale_fill_manual(name = "Temperature",
                      values = palette,
                      labels = c("18.5°C", 
                                 "23.5°C", 
                                 "28.5°C", 
                                 "33.5°C", 
                                 "38.5°C")) +
    scale_color_manual(values = palette) + 
    labs(x = "Ethanol (mM)",
         y = expression("Reaction Velocity (Abs/sec)")) + 
    theme_classic()
} 
# End function -----------------------------------------------------------------


# ------------------------------------------------------------------------------
# Function: plot_km
# Description: Plot km values with std.error bars
# Inputs: tibble with Km estimate and error
# Outputs: plot

require(tidyverse)

plot_km <- function(data, 
                    palette = brewer.pal(4, "Set2"),
                    legend.labels = c( "VT-ECK 8", "Chiapas", "Fast", "Slow")) {
  ggplot(data, aes(x = temp_group, y = estimate_Km, fill = enzyme)) +
    geom_errorbar(aes(ymin = estimate_Km - std.error_Km, 
                      ymax = estimate_Km + std.error_Km), 
                  width = .2,
                  position = position_dodge(0.5)) +
    geom_point(pch = 21, position = position_dodge(0.5), size = 6) +
    labs(x = "Temperature (°C)",
         y = expression("K"[m]*"  (mM EtOH)")) +
    scale_fill_manual(name = "ADH enzyme",
                      values = palette,
                      labels = legend.labels) +
    theme_classic()
} 
# End function -----------------------------------------------------------------