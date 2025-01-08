#BB
#10/30/24
#Trying some regressions to begin with

library(ggplot2)
library(mvabund)
library(readxl)
library(dplyr)
library(tidyr)

setwd("C:/Users/beabo/OneDrive/Documents/NAU/Classes Archived/Australia Co-Invasions Course 2023/Aus-Invasions-2023-Course/Merged_data")

ds <- read.csv("Aus23_allData_19Nov24.csv")

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


hist(ds$ergosterol) #Looks more like a gamma dist...

mod <- glm(ergosterol ~ ECM_abund_soil+TreeSpecies, family = "Gamma", ds)
summary(mod)

lin_mod <- lm(ergosterol~ ECM_abund_soil+TreeSpecies, ds)
summary(lin_mod)

#Interaction term is not significant. ECM abundance in soil seems to be a part of the ergosterol. Ergosterol higher around pines. 

plot <- ds %>%
  ggplot(aes(y = ergosterol, x = ECM_abund_soil, color = TreeSpecies)) +
  geom_point() +
  scale_color_manual(values = c("#EAC435", "#345995"))+
  geom_smooth(aes(group = TreeSpecies), method = "glm", method.args = list(family = "Gamma"),) +
  labs(y = "Ergosterol", x  = "ECM Abundance (soil)")

ggsave("C:/Users/beabo/OneDrive/Documents/NAU/Classes Archived/Australia Co-Invasions Course 2023/Aus-Invasions-2023-Course/Plots/Analyses/ergosterol_vs_ecm_tree.png", plot = plot, dpi = 500, width = 5, heigh = 3, units = "in")


#Do some sort of AIC-based model selection with all predictors
#Response variables; Total CNP, decomp, 
#Decomp 
#Dredge package


#Seems like fungi in roots respond more to treespecies than bacteria do, while bacteria in soil respond more to treespecies than fungi do. How to evaluate that?

#Root, soil
#Bac vs. fung

long_data <- ds %>%
  pivot_longer(
    cols = ends_with("_soil") | ends_with("_root"),  # Select columns to pivot
    names_to = c("Microbe_Type", "Metric", "Sample_Location"),  # New columns for Microbe_Type and Sample_Location
    names_pattern = "([A-Za-z]+)_(.*)_(soil|root)",  # Split the name into Microbe_Type and Sample_Location
    values_to = "Value"  # Column for the values
  )  %>%
  mutate(
    Bac_or_Fung = case_when(
      Microbe_Type %in% c("Bac", "OtherBacteria") ~ "Bacteria",
      Microbe_Type %in% c("AMF", "ECM", "Fungal", "Endo", "Sap") ~ "Fungi",
      TRUE ~ Microbe_Type  # Optional: set to NA for other values
    )
  ) %>%
  pivot_wider(
    names_from = Metric,  # Keep metrics as separate columns
    values_from = Value,  # Fill new columns with corresponding values
    values_fill = list(Value = NA)  # Fill NAs for missing values
  )
  #relocate(colnames(long_data)[51:ncol(long_data)])
  
long_data %>%
  group_by(Bac_or_Fung, Sample_Location, TreeSpecies)%>%
  summarize(n = n())
  
View(long_data%>%
  filter(Bac_or_Fung == "Fungi" & Sample_Location == "soil" & TreeSpecies == "Pine")) #not sure what is up with this

mod <- lm(richness ~ Bac_or_Fung:Sample_Location:TreeSpecies, long_data)
summary(mod)
#Bacteria in root of euc: +richness
#fungi in root of euc: -richness
#Bac in soil around euc: +richness, more so than in root (x2)
#Fungi in soil around euc: +richness (flipped sign from root)

#Bac in root of pine: +richness
#Fungi in root of pine: -richness, bigger response than in euc roots
#Bac in soil around pine: + richness (more so than in root, 1.5x)

long_data %>%
  filter(Bac_or_Fung %in% c("Bacteria", "Fungi"))%>%
  ggplot(aes(y = Richness, x = TreeSpecies, fill =Sample_Location ))+
  geom_boxplot()+
  facet_grid(~Bac_or_Fung)
#not much to visualize

mod <- lm(abund ~ Bac_or_Fung:Sample_Location + TreeSpecies, long_data)
summary(mod)

mod <- lm(Richness ~ Bac_or_Fung*Sample_Location*TreeSpecies + soil_moisture, long_data)
