#BSA

# Load packages
require(tidyverse)

# Bradford Assay -----

# Load data
df <- read_csv(here::here("ox_stress/sod/data/pilot_bsa_060221.csv"))

# Standard data
df_std <- df %>%
  filter(Sample == "Standard")

# Plot Standard Curve
df_std %>%
  ggplot(aes(x = Conc, y = Abs)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "grey50") +
  theme_classic()

# Run Linear Regression on data 
fit <- lm(df_std, formula = Abs ~ Conc)

summary(fit)$r.squared
y_int <- fit$coefficients[1]
slope <- fit$coefficients[2]


# Calculate dilutions 
final_protein_conc <- 1000
final_soln_vol <- 50

df_calc <- df %>%
  filter(Sample !="Standard") %>%
  mutate(Conc_calc = (Abs-y_int)/slope) %>%
  group_by(Sample) %>%
  summarise(Conc_calc = mean(Conc_calc)) %>%
  mutate(vol_soln = (final_protein_conc*final_soln_vol)/Conc_calc,
         vol_buffer = final_soln_vol - vol_soln)

# SOD Activity Assay ----

# Load data
df <- read_csv(here::here("ox_stress/sod/data/pilot_sod_activity_060221.csv"))

# Standard Curve
df_std <- df %>%
  filter(Sample == "Standard") %>%
  group_by(Sample, Activity) %>%
  summarise(Abs = mean(Abs)) %>%
  mutate(Ratio = Abs[Activity == 0]/Abs)

# Run Linear Regression
fit <- lm(df_std, formula = Ratio ~ Activity)

summary(fit)$r.squared

y_int <- fit$coefficients[1]
slope <- fit$coefficients[2]


# Plot Standard Curve
df_std %>%
  ggplot(aes(x = Activity, y = Ratio)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "grey50") +
  theme_classic()

# Calculate Activity
df_calc <- df %>%
  filter(!is.na(Total_protein_conc) | Sample == "Standard" & Activity == 0) %>%
  group_by(Sample, Total_protein_conc) %>%
  #summarise(Abs = mean(Abs)) %>%
  ungroup() %>%
  mutate(Ratio = Abs[Sample == "Standard"]/Abs,
         sod_activity_U_mL = ((Ratio - y_int)/slope) * (0.230/0.01))

df_calc %>%
  filter(Sample != "Standard") %>%
  ggplot() +
  geom_point(aes(x = Total_protein_conc, y = Ratio, color = Sample)) + 
  theme_classic()
