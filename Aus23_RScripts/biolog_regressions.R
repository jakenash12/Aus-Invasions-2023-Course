#BB
#10/30/24
#Trying some regressions to begin with

library(ggplot2)
library(mvabund)
library(readxl)
library(dplyr)

ds <- read_xlsx("Aus-Invasions-2023-Course/Merged_data/Aus23_master_pooled.xlsx")

theme_set(theme_bw())

# Litter decomp -----------------------------------------------------------
#How does ECM affect litter decomp in the two types of soils?

mod <- lm(litter_depth ~ soil_moisture + ECM_richness_soil*TreeSpecies, ds)
summary(mod)
#Litter depth is higher when tree species is pine. Litter depth higher when soil moisture is higher.  Interaction term of ECM richness of soil and tree species pine leads causes small decreases in litter depth.

plot <- ds %>%
  ggplot(aes(y = litter_depth, x = ECM_richness_soil)) +
  geom_point(aes(color = soil_moisture)) +
  geom_smooth(aes(group = TreeSpecies), method = "glm", color = "black") +
  scale_color_gradient(low = "lightblue", high = "darkblue", name = "Soil Moisture")+
  facet_grid(~TreeSpecies)+
  labs(y = "Litter Depth (cm)", x  = "ECM Richness (soil)")

ggsave("Aus-Invasions-2023-Course/Plots/Analyses/litterdepth_vs_ecm_moisture_tree.png", plot = plot, dpi = 500, width = 5, heigh = 3, units = "in")

mod <- lm(ergosterol ~ ECM_abund_soil + TreeSpecies, ds)
summary(mod)
#Interaction term is not significant. ECM abundance in soil seems to be a massive part of the ergosterol. Ergosterol also higher around pines. 

plot <- ds %>%
  ggplot(aes(y = ergosterol, x = ECM_abund_soil, color = TreeSpecies)) +
  geom_point() +
  scale_color_manual(values = c("#EAC435", "#345995"))+
  geom_smooth(aes(group = TreeSpecies), method = "glm") +
  labs(y = "Ergosterol", x  = "ECM Abundance (soil)")

ggsave("Aus-Invasions-2023-Course/Plots/Analyses/ergosterol_vs_ecm_tree.png", plot = plot, dpi = 500, width = 5, heigh = 3, units = "in")


#Do some sort of AIC-based model selection with all predictors
#Response variables; Total CNP, decomp, 
#Decomp 
#Dredge package

global.model <- 