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
                 names_to = "conc_well", 
                 values_to = "abs") %>%
    mutate(conc = as.numeric(str_replace_all(conc_well, "_.", ""))) %>%
    mutate(file = str_remove(file, ".csv")) %>%
    separate(file, into = c("enzyme", "temp_group", "date"), sep = "_")
  
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
    group_by(enzyme, conc_well, temp_group) %>%
    nest() %>%
    mutate(
      fit = map(data, ~ lm(abs ~ time, data = .x)),
      tidied = map(fit, tidy),
      glanced = map(fit, glance)
    ) %>% 
    unnest(tidied) %>%
    unnest(glanced, .name_repair = "minimal") %>%
    filter(term == "time") %>%
    select(enzyme, temp_group, conc_well, estimate, r.squared) %>%
    rename(velocity = estimate) %>%
    mutate(conc = as.numeric(str_replace_all(conc_well, "_.", ""))) %>%
    filter(conc != 0) %>%
    arrange(conc) # %>%
    # group_by(temp_group, conc) %>%
    # summarize(velocity = mean(velocity))
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
    group_by(enzyme, temp_group) %>%
    nest() %>%
    mutate(
      fit = map(data, ~ nls(formula = 
                              velocity ~ SSmicmen(conc, Vmax, Km), data = .), 
                data = .x),
      tidied = map(fit, tidy)
    ) %>%
    unnest(tidied) %>%
    mutate(r.squared = map(fit, Rsq)) %>%
    unnest(r.squared)
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
                   color = as.factor(conc_well))) +
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap(~ temp_group) +
    labs(x = "Time (s)", 
         y = "Absorbance @ 340 nm", 
         title = str_to_title(unique(data$enzyme))) +
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
                          palette = c("skyblue", "blue", "purple", 
                                      "red", "darkred"),
                          legend_labels = c("18.5°C", "23.5°C", "28.5°C", 
                                            "33.5°C", "38.5°C")) {
  
  ggplot(data, aes(x = conc, y = velocity, fill = temp_group)) +
    geom_point(pch = 21, size = 3, color = "grey50", alpha = 0.75) +
    geom_smooth(method = "nls", 
                formula = y ~ SSmicmen(x, Vmax, Km),
                data = data,
                se = FALSE,
                aes(color = temp_group), 
                show.legend = FALSE) +
    # scale_fill_manual(name = "Temperature",
    #                   values = palette,
    #                   labels = legend_labels) +
    # scale_color_manual(values = palette) + 
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
                    legend_labels = c( "VT-ECK 8", "Chiapas", "Fast", "Slow")) {
  ggplot(data, aes(x = temp_group, y = estimate_Km, fill = enzyme)) +
    geom_errorbar(aes(ymin = estimate_Km - std.error_Km, 
                      ymax = estimate_Km + std.error_Km), 
                  width = .2,
                  position = position_dodge(0.5)) +
    geom_point(pch = 21, position = position_dodge(0.5), size = 6) +
    labs(x = "Temperature (°C)",
         y = expression("K"[m]*"  (mM EtOH)")) +
    # scale_fill_manual(name = "ADH enzyme",
    #                   values = palette,
    #                   labels = legend_labels) +
    theme_classic()
} 
# End function -----------------------------------------------------------------



# ------------------------------------------------------------------------------
# Function: Rsq
# Description: Function to calculate the multiple R-squared and the adjusted 
#              R-squared from a fitted model via lm or aov, i.e., linear models. 
#              For a model fitted via nls, nonlinear models, the pseudo 
#              R-squared is returned. Rsq function from soilphysics package
# Inputs: a fitted model via lm or aov or via nls, nonlinear models
# Outputs: output_description

Rsq <- function(model) {
  if (!inherits(model, c("lm", "aov", "nls")))
    stop ("'Rsq' is only applied to the classes: 'lm', 'aov' or 'nls'.")
  if (inherits(model, c("glm", "manova", "maov", "mlm")))
    stop("'Rsq' is not applied to an object of this class!")
  
  pred <- predict(model)
  n <- length(pred)
  res <- resid(model)
  w <- weights(model)
  if (is.null(w)) w <- rep(1, n)
  rss <- sum(w * res ^ 2)
  resp <- pred + res
  center <- weighted.mean(resp, w)
  if (inherits(model, c("lm", "aov"))) {
    r.df <- model$df.residual
    int.df <- attr(model$terms, "intercept")
    if (int.df) {
      mss <- sum(w * scale(pred, scale = FALSE)^2)
    } else {
      mss <- sum(w * scale(pred, center = FALSE,
                           scale = FALSE)^2)
    }
    r.sq <- mss / (mss + rss)
    adj.r.sq <- 1 - (1 - r.sq) * (n - int.df) / r.df
    out <- list(R.squared = r.sq, adj.R.squared = adj.r.sq)
  } else {
    r.df <- summary(model)$df[2]
    int.df <- 1
    tss <- sum(w * (resp - center)^2)
    r.sq <- 1 - rss/tss
    adj.r.sq <- 1 - (1 - r.sq) * (n - int.df) / r.df
    out <- list(pseudo.R.squared = r.sq,
                adj.R.squared = adj.r.sq)
  }
  class(out) <- "Rsq"
  return(out$pseudo.R.squared)
}
# End function -----------------------------------------------------------------

# ------------------------------------------------------------------------------
# Function: read_tidy_nadh_csv
# Description: Reads a .csv with the nadh std curve values and tidys the data
# Inputs: file -- .csv file location
# Outputs: tidy tibble with NADH, mins, and abs columns 

require(tidyverse)

read_tidy_nadh_csv <- function(file) {
  df <- read_csv(file)
  
  tidy_df <- df %>%
    pivot_longer(cols = contains("min"), names_to = "mins", values_to = "abs") %>%
    mutate(mins = as.numeric(str_remove_all(mins, "_min"))) 
  
  return(tidy_df)
  
} 
# End function -----------------------------------------------------------------


# ------------------------------------------------------------------------------
# Function: read_tidy_adh_abs_csv
# Description: Reads and tidys ADH absorbance data
# Inputs: file -- .csv location
# Outputs: tidy tibble with EtOH mM, sample_vol, n_flies, genotype, mins and abs
#            columns

require(tidyverse)

read_tidy_adh_abs_csv <- function(file) {
  df <- read_csv(file)
  
  df <- df %>%
    pivot_longer(cols = contains("min"), names_to = "mins", values_to = "abs") %>%
    mutate(mins = as.numeric(str_remove_all(mins, "_min")))
  
  return(df)
} 
# End function -----------------------------------------------------------------


# ------------------------------------------------------------------------------
# Function: plot_std_curve
# Description: Plots the NADH std curve
# Inputs: Tidy tibble with NADH and abs cols -- abs average or one time point
# Outputs: plot with linear regression curve

plot_std_curve <- function(dat) {
  ggplot(data = dat, aes(x = NADH, y = abs)) +
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
} 
# End function -----------------------------------------------------------------


# ------------------------------------------------------------------------------
# Function: std_curve
# Description: Get intercept and slope values of std curve for downstream calcs
# Inputs: data.frame with NADH and abs columns
# Outputs: output_description

std_curve <- function(dat) {
  int <- as.numeric(lm(abs ~ NADH, dat)$coefficients[1])
  slope <- as.numeric(lm(abs ~ NADH, dat)$coefficients[2])
  
  return(c(int = int, slope = slope))
} 
# End function -----------------------------------------------------------------


# ------------------------------------------------------------------------------
# Function: calc_adh_activity
# Description: description
# Inputs: input_description
# Outputs: output_description

calc_adh_activity <- function(dat, std_curve_lm, t_0 = 0) {
  dat %>%
    mutate(nmol_NADH = (abs - std_curve_lm["int"]) / std_curve_lm["slope"]) %>%
    mutate(deltaAbs = abs - abs[mins == t_0],
           deltaNADH = nmol_NADH - nmol_NADH[mins == t_0],
           deltaTime = mins - t_0) %>%
    filter(deltaTime != 0) %>%
    mutate(dilution_factor = 150/sample_vol) %>%
    mutate(mU_mL = (deltaNADH/(deltaTime*(sample_vol/1000)))*dilution_factor) %>%
    ungroup(EtOH) %>%
    mutate(mU_mL = mU_mL - mU_mL[EtOH == 0])
} 
# End function -----------------------------------------------------------------


# ------------------------------------------------------------------------------
# Function: plot_sample_abs_time
# Description: description
# Inputs: input_description
# Outputs: output_description

plot_sample_abs_time <- function(dat) {
  ggplot(data = dat, 
         aes(x = mins, y = abs, 
             color = EtOH)) +
    geom_smooth(method='lm', formula = y ~ x, se = FALSE) +
    geom_point() +
    theme_classic() +
    labs(y = "Absorbance @ 450 nm", x = "Time (mins)") 
} 
# End function -----------------------------------------------------------------


# ------------------------------------------------------------------------------
# Function: plot_mm_curve_biovision
# Description: Plot the Michaelis-Menten curves
# Inputs: tibble with conc and velocity columns
# Outputs: plot

require(tidyverse)

plot_mm_curve_biovision <- function(data) {
  
  ggplot(data, aes(x = as.numeric(as.character(EtOH)), y = mU_mL)) +
    geom_point(pch = 21, size = 3, color = "grey50", alpha = 1) +
    geom_smooth(method = "nls", 
                formula = y ~ SSmicmen(x, Vmax, Km),
                data = data,
                se = FALSE,
                show.legend = FALSE) +
    labs(x = "Ethanol (mM)",
         y = expression("ADH activity (mU/mL)")) + 
    theme_classic()
} 
# End function -----------------------------------------------------------------




# ------------------------------------------------------------------------------
# Function: plot_mm_curve_biovision
# Description: Plot the Lineweaver-Burke
# Inputs: tibble with conc and velocity columns
# Outputs: plot

require(tidyverse)

plot_mm_curve_lb <- function(data) {
  
  ggplot(data, aes(x = 1/as.numeric(as.character(EtOH)), y = 1/mU_mL)) +
    geom_point(pch = 21, size = 3, color = "grey50", alpha = 1) +
    geom_smooth(method = "lm", 
                formula = y ~ x,
                data = data,
                se = FALSE,
                show.legend = FALSE) +
    labs(x = "1/[Ethanol (mM)]",
         y = expression("1/(ADH activity (mU/mL))")) + 
    theme_classic()
} 
# End function -----------------------------------------------------------------





# ------------------------------------------------------------------------------
# Function: calc_Km_Vm_biovision
# Description: Calculate Vmax and Km for Michaelis–Menten curve
# Inputs: tibble with velocity estimates for each conc and temp_group
# Outputs: tibble with estimates for Vmax and Km

require(tidyverse)
require(broom)

calc_Km_Vm_biovision <- function(data) {
  data %>%
    mutate(EtOH = as.numeric(as.character(EtOH))) %>%
    filter(EtOH != 0) %>%
    group_by(genotype, temp) %>%
    nest() %>%
    mutate(
      fit = map(data, ~ nls(formula = 
                              mU_mL ~ SSmicmen(EtOH, Vmax, Km), data = .), 
                data = .x),
      tidied = map(fit, tidy)
    ) %>%
    unnest(tidied) %>%
    mutate(r.squared = map(fit, Rsq)) %>%
    unnest(r.squared)
} 
# End function -----------------------------------------------------------------