# ------------------------------------------------------------------------------
# Survival versus gene expression for oxidative stress genes
# April 12, 2021
# TS O'Leary
# ------------------------------------------------------------------------------

# Load packages
require(tidyverse)
require(DESeq2)

# Load phenotype data
lt_df <- read_delim(here::here("ox_stress/lockwood_2018_LT50_data.txt"), 
                    delim = "\t")

survival <- read_csv(here::here("ox_stress/lockwood_2018_raw_survival.csv"))

# Load expression data
ddslrtreg <- readRDS(here::here("ox_stress/transcript_shiny/ddslrtreg.rds"))


# Quick survival plot
survival %>% 
  group_by(temperature, genotype) %>%
  summarise(mean_survival = mean(survival)) %>%
  ggplot() +
  geom_point(aes(x = temperature, y = mean_survival, color = genotype))

# Summarize mean & filter to 
df <- survival %>% 
  group_by(temperature, genotype, region) %>%
  summarise(mean_survival = mean(survival)) %>%
  filter(region == "tropical" | str_detect(genotype, "VTECK_")) %>%
  mutate(genotype = str_remove_all(genotype, "ECK_"),
         genotype = str_replace_all(genotype, "MU", "BO"),
         genotype = str_replace_all(genotype, "CH", "CP")) %>%
  arrange(genotype, temperature)


# Correlate LT50 and expression data -------------------------------------------

transcript_id <- "FBtr0304051"

# Plot each
d <- plotCounts(ddslrtreg, 
                gene = transcript_id, 
                intgroup = c("pop", "temp"), 
                returnData = TRUE) %>%
  left_join(lt_df, by = c("pop" = "Pop"))


# Plot simple expression v LT50 plot 
x %>%
  ggplot(aes(x = Embryo_LT50, y = count)) +
  geom_boxplot(aes(fill = Pop), color = "grey50", alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = transcript_id) +
  theme_classic() +
  facet_wrap(~ temp)


