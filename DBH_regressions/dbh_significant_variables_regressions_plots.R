library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)


Aus23_allData_19Nov24 <- read.csv("Merged_data/Aus23_allData_19Nov24.csv") %>% select(-X)
variables.signif.take2<-c("Tween_40_BiologDay5", "L.Serine_BiologDay5", "Glycyl.L.Glutamic_Acid_BiologDay5",
                          "i.Erythritol_BiologDay5", "Putrescine_BiologDay5", "L.Asparagine_BiologDay5", 
                          "pine_litter_prop", "Euc_leafLitter_percP", "Euc_leafLitter_totP",              
                          "Litter_OLayer_ergosterol", "Litter_avg_ergosterol","euc_litter_prop", "Path_abund_root", "ergosterol","soil_moisture",
                          "Bac_Shannon_soil")



# Filter and pivot data to long format
data_long_signif <- Aus23_allData_19Nov24 %>%
  select(all_of(c("dbh_cm", "TreeSpecies", variables.signif.take2))) %>%
  pivot_longer(
    cols = -c(dbh_cm, TreeSpecies), 
    names_to = "Response", 
    values_to = "Value"
  )

# Plot regressions for each variable against dbh_cm
ggplot(data_long_signif, aes(x = dbh_cm, y = Value, color = TreeSpecies)) +
  geom_point(size = 2, alpha = 0.7) +  # Scatterplot points
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line
 # geom_smooth(method = "loess", se = TRUE, linetype = "solid", color = "black") +  # LOESS trend line with confidence intervals
  facet_wrap(~ Response, scales = "free_y") +  # Facet by response variable
  labs(title = "Regression of Significant Variables Against DBH",
       x = "DBH (cm)",
       y = "Response Variable Value") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 10, face = "bold"),  # Facet label formatting
    axis.title = element_text(size = 12)
  )

ggsave("Regression_Plots_DBH_vs_Signif_Variables_Faceted.png", width = 14, height = 8, dpi = 300)



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


