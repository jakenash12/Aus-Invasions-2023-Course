# Regressions of variables that differ between Euc vs Pine by tree DBH
# Graphs started Dec 2024 by Caitlín Dagg
# Stats Jan 2025 by Corinne V

library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(ggpmisc)
library(car)
library(lmerTest)

setwd("Aus-Invasions-2023-Course/")

fulldat <- read.csv("Merged_data/Aus23_allData_18Jan25.csv", row.names=1)

################# RUN REGRESSIONS ########################################################

# select out the variables we want to run regressions on
vars <- fulldat[c(4:77,80:84,150:206)]

# write a function to run the regressions:
run_regressions <- function(variable_name, data) {
  formula <- as.formula(paste(variable_name, "~ dbh_cm*TreeSpecies"))
  mod <- lm(formula, data = data)
  mod.summ <- summary(mod)
  mod.R2 <- mod.summ$r.squared
  mod.anova <- car::Anova(mod, test.statistic = "F")
  c(dbh_p = mod.anova[1,4],
    tree_p = mod.anova[2,4],
    int_p = mod.anova[3,4],
    dbh_f = mod.anova[1,3],
    tree_f = mod.anova[2,3],
    int_f = mod.anova[3,3],
    R2 = mod.R2)
}

# apply this function to all variables
variables <- setdiff(names(vars), c("dbh_cm", "TreeSpecies"))
lm.results <- lapply(variables, run_regressions, data = vars)

# convert into dataframe
lm.pvals <- do.call(rbind, lm.results)
lm.pvals <- as.data.frame(lm.pvals)
lm.pvals$variable <- variables

# pull out all variables where dbh or the interaction are significant
lm.signif <- subset(lm.pvals, dbh_p < 0.05 | int_p < 0.05)

# select variables where the interaction is significant OR both tree and DBH are significant (excluding redundnat variables like diversity metrics)
# AMF_abund_soil has a significant interaction but it's driven entirely by one point - not including
select.vars <- c("dbh_cm", "TreeSpecies", "Bac_Shannon_soil", "ECM_richness_soil", "AMF_abund_root", "Endo_abund_root",
                 "soil_moisture", "pine_litter_prop", "euc_litter_prop", "Tween_40_BiologDay5", "a.Cyclodextrin_BiologDay5")
select.dat <- fulldat[select.vars]

####### Check linear model assumptions ##########
mod <- lm(pine_litter_prop ~ dbh_cm*TreeSpecies, data = fulldat) # change response variable to what you want to look at
summary(mod)
anova(mod)
plot(mod) # plots of residuals and qqnorm 


############################# PLOTS #############################

# set order of plots and variable names
var.order <- c("Tween_40_BiologDay5", "a.Cyclodextrin_BiologDay5", "Bac_Shannon_soil", "ECM_richness_soil", "AMF_abund_root", "Endo_abund_root",
               "soil_moisture", "pine_litter_prop", "euc_litter_prop")
y_labels <- c("Tween 40 utilization", "Cyclodextrin utilization", "Soil bacterial Shannon diversity", "Soil EMF richness", 
              "Root AMF relative abundance", "Root endophytic fungi relative abundance",
              "Soil moisture (proportion)", "Proportion pine litter", "Proportion eucalypt litter")


# Filter and pivot data to long format
data_long_signif <- select.dat %>%
  pivot_longer(
    cols = -c(dbh_cm, TreeSpecies), 
    names_to = "Response", 
    values_to = "Value")
data_long_signif$Response <- factor(data_long_signif$Response, levels = var.order)

# Plot regressions for each variable against dbh_cm
ggplot(data_long_signif, aes(x = dbh_cm, y = Value, color = factor(TreeSpecies, c("Pine", "Eucalyptus")))) +
  geom_point(size = 2, alpha = 0.7) +  # Scatterplot points
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, aes(linetype = TreeSpecies)) +  # Linear regression line
  facet_wrap(~ Response, scales = "free_y", labeller = as_labeller(y_labels), strip.position = "left") +
  labs(x = "",
       y = "") +
  scale_color_manual(values = c("#000DCC", "#FF9900")) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.spacing.x = unit(4, "lines"),
    panel.spacing.y = unit(1, "lines"),
    strip.text = element_blank(),  # Facet label formatting
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.title = element_text(size = 12),
    axis.text = element_text(colour="black", size = 12),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()
  )

ggsave("DBH_regressions/DBH_regressions_formatted.png", width = 12.5, height = 8, dpi = 300)







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

understory.lm <- lm(Understory_richness_perTree ~ dbh_cm*TreeSpecies, data = fulldat)
Anova(understory.lm, test.statistic = "F")
plot(understory.lm)

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



