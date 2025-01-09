# Regressions of variables that differ between Euc vs Pine by tree DBH
# Graphs started Dec 2024 by Caitlín Dagg
# Stats Jan 2025 by Corinne V


library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(ggpmisc)
library(car)
install.packages("ggpmisc")


Aus23_allData_19Nov24 <- read.csv("Aus-Invasions-2023-Course/Merged_data/Aus23_allData_19Nov24.csv") %>% select(-X)
variables.signif.take2<-c("Bac_Shannon_soil","Path_abund_root", "ergosterol", "Tween_40_BiologDay5", "L.Serine_BiologDay5", 
                          "Glycyl.L.Glutamic_Acid_BiologDay5",
                          "i.Erythritol_BiologDay5", "Putrescine_BiologDay5", "L.Asparagine_BiologDay5", "soil_moisture",
                          "pine_litter_prop", "euc_litter_prop", "Euc_leafLitter_percP",           
                          "Litter_OLayer_ergosterol", "Litter_avg_ergosterol")

############################# PLOTS ############################# 

# Filter and pivot data to long format
data_long_signif <- Aus23_allData_19Nov24 %>%
  select(all_of(c("dbh_cm", "TreeSpecies", variables.signif.take2))) %>%
  pivot_longer(
    cols = -c(dbh_cm, TreeSpecies), 
    names_to = "Response", 
    values_to = "Value")
data_long_signif$Response <- factor(data_long_signif$Response, levels = variables.signif.take2)

# Plot regressions for each variable against dbh_cm
ggplot(data_long_signif, aes(x = dbh_cm, y = Value, color = factor(TreeSpecies, c("Pine", "Eucalyptus")))) +
  geom_point(size = 2, alpha = 0.7) +  # Scatterplot points
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, aes(linetype = TreeSpecies)) +  # Linear regression line
  facet_wrap(~ Response, scales = "free_y") +  # Facet by response variable
  labs(x = "DBH (cm)",
       y = "Response Variable Value") +
  scale_color_manual(values = c("#000DCC", "#FF9900")) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.spacing.x = unit(2, "lines"),
    strip.text = element_text(size = 10, face = "bold"),  # Facet label formatting
    axis.title = element_text(size = 12),
    axis.text = element_text(colour="black"),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()
  )

ggsave("Aus-Invasions-2023-Course/DBH_regressions/Regression_Plots_DBH_vs_Signif_Variables_Faceted_Stats_v2.png", width = 14, height = 8, dpi = 300)


# Loop through each significant variable and plot separately
for (var in variables.signif.take2) {
  
  # Subset data for the current variable
  data_subset <- data_long_signif %>% filter(Response == var)
  
  # Create the plot
  p <- ggplot(data_subset, aes(x = dbh_cm, y = Value, color = TreeSpecies)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
 #   geom_smooth(method = "loess", se = TRUE, linetype = "solid") +
    labs(title = paste("Regression of", var, "Against DBH"),
         x = "DBH (cm)",
         y = var) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12)
    )
  
  # Save the plot as a PNG file
  ggsave(paste0( "DBH_vs_", var, ".png"), plot = p, width = 8, height = 6, dpi = 300)
}

############################# STATS ############################# 
hist(Aus23_allData_19Nov24$dbh_cm)

tween.lm <- lm(Tween_40_BiologDay5 ~ dbh_cm*TreeSpecies, data = Aus23_allData_19Nov24)
Anova(tween.lm, test.statistic = "F") 
plot(tween.lm) # dbh, tree

moisture.lm <- lm(soil_moisture ~ dbh_cm*TreeSpecies, data = Aus23_allData_19Nov24)
Anova(moisture.lm, test.statistic = "F") 
plot(moisture.lm) # dbh, tree, int

Lserine.lm <- lm(L.Serine_BiologDay5 ~ dbh_cm*TreeSpecies, data = Aus23_allData_19Nov24)
Anova(Lserine.lm, test.statistic = "F")
plot(Lserine.lm) # tree marginal

glycyl.lm <- lm(Glycyl.L.Glutamic_Acid_BiologDay5 ~ dbh_cm*TreeSpecies, data = Aus23_allData_19Nov24)
Anova(glycyl.lm, test.statistic = "F")
plot(glycyl.lm) # tree
# residuals are eh...

iEry.lm <- lm(i.Erythritol_BiologDay5 ~ dbh_cm*TreeSpecies, data = Aus23_allData_19Nov24)
Anova(iEry.lm, test.statistic = "F")
plot(iEry.lm) # tree

putre.lm <- lm(Putrescine_BiologDay5 ~ dbh_cm*TreeSpecies, data = Aus23_allData_19Nov24)
Anova(putre.lm, test.statistic = "F")
plot(putre.lm) # tree

LAsp.lm <- lm(L.Asparagine_BiologDay5 ~ dbh_cm*TreeSpecies, data = Aus23_allData_19Nov24)
Anova(LAsp.lm, test.statistic = "F")
plot(LAsp.lm) # tree

pinelit.lm <- lm(pine_litter_prop ~ dbh_cm*TreeSpecies, data = Aus23_allData_19Nov24)
Anova(pinelit.lm, test.statistic = "F")
plot(pinelit.lm) # tree, dbh

eucP.lm <- lm(Euc_leafLitter_percP ~ dbh_cm*TreeSpecies, data = Aus23_allData_19Nov24)
Anova(eucP.lm, test.statistic = "F")
plot(eucP.lm) # tree

Oergosterol.lm <- lm(Litter_OLayer_ergosterol ~ dbh_cm*TreeSpecies, data = Aus23_allData_19Nov24)
Anova(Oergosterol.lm, test.statistic = "F")
plot(Oergosterol.lm) # tree

avgLitErg.lm <- lm(Litter_avg_ergosterol ~ dbh_cm*TreeSpecies, data = Aus23_allData_19Nov24)
Anova(avgLitErg.lm, test.statistic = "F")
plot(avgLitErg.lm) # tree

euclit.lm <- lm(euc_litter_prop ~ dbh_cm*TreeSpecies, data = Aus23_allData_19Nov24)
Anova(euclit.lm, test.statistic = "F")
plot(euclit.lm) # dbh, tree

rootpath.lm <- lm(Path_abund_root ~ dbh_cm*TreeSpecies, data = Aus23_allData_19Nov24)
Anova(rootpath.lm, test.statistic = "F")
plot(rootpath.lm) # tree
# residuals iffy

soilErg.lm <- lm(ergosterol ~ dbh_cm*TreeSpecies, data = Aus23_allData_19Nov24)
Anova(soilErg.lm, test.statistic = "F")
plot(soilErg.lm) # tree

bacshannon.lm <- lm(Bac_Shannon_soil ~ dbh_cm*TreeSpecies, data = Aus23_allData_19Nov24)
Anova(bacshannon.lm, test.statistic = "F")
plot(bacshannon.lm) # dbh, tree




# specify y axis labels for each variable
Ylabels <- c("Tween40 metabolism", "L-Serine metabolism", "Glycyl-L-Glutamic Acid metabolism",
             "i-Erythritol metabolism", "Putrescine metabolism", "L-Asparagine metabolism", 
             "Proportion pine litter", "Proportion eucalypt litter", "Eucalpyt leaf litter %P",           
             "Litter O layer ergosterol", "All litter layer ergosterol", "Soil ergosterol content", "Root fungal pathogen relative abundance", "Soil moisture",
             "Soil bacterial Shannon diversity")

variables.signif.take2<-c("Tween_40_BiologDay5", "L.Serine_BiologDay5", "Glycyl.L.Glutamic_Acid_BiologDay5",
                          "i.Erythritol_BiologDay5", "Putrescine_BiologDay5", "L.Asparagine_BiologDay5", 
                          "pine_litter_prop", "Euc_leafLitter_percP", "Euc_leafLitter_totP",              
                          "Litter_OLayer_ergosterol", "Litter_avg_ergosterol","euc_litter_prop", "Path_abund_root", 
                          "ergosterol","soil_moisture",
                          "Bac_Shannon_soil")



