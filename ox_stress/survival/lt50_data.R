# ------------------------------------------------------------------------------
# LT50 of Natural Lines for redox experiments
# April 12, 2021
# TS O'Leary
# ------------------------------------------------------------------------------

# Load packages
require(tidyverse)

# Load data
lt_df <- read_delim(here::here("ox_stress/lockwood_2018_LT50_data.txt"), 
                    delim = "\t")
  

# Embryo LT50 Phenotype data visualize
lt_df %>%
  filter(Region == "Tropical" | StateCountry == "Vermont_USA") %>%
  ggplot(aes(x = Region, 
             y = Embryo_LT50, 
             fill = Region, 
             label = Pop)) +
  geom_boxplot() +
  geom_jitter(color = "black", fill = "grey50",
              shape = 21,
              size = 3,
              alpha = 0.85,
              position = position_jitter(0)) +
  ggrepel::geom_label_repel(alpha = 0.7) +
  scale_fill_manual(values = c("#DDCC77", "#CC6677")) +
  ylab(expression("LT"[50]*" (°C)")) +
  labs(title = "Embryonic thermal tolerance") +
  theme_classic(base_size = 18) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))

# Min Vermont Temperate
lt_df %>%
  filter(StateCountry == "Vermont_USA") %>%
  group_by(Region) %>%
  top_n(-4, Embryo_LT50)

# Max Tropical
lt_df %>%
  filter(Region == "Tropical") %>%
  group_by(Region) %>%
  top_n(2, Embryo_LT50)

