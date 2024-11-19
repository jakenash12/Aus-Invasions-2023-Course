#AU Manuscript Analyses
#Generating a forest plot to compare effect sizes of each response variable in Pine vs Euc
#This is the new forest plot that accounts for background correction in Biolog Ecoplates & gets rid of redundant/unwanted response variables 

setwd("/Users/elenaleander/Desktop/AU_manuscript/") 
library(metafor)
library(dplyr)
library(tidyr)

#update masterdata with Biolog Ecoplate background correction (ie subtract water values from substrate values) 
AU_masterdata <- read.csv("/Users/elenaleander/Desktop/AU_manuscript/Aus23_allData_19Nov24.csv", header = TRUE)

########Cohen's d, without response variable grouping#########
#Use when sample sizes are relatively equal (Euc & Pine both = 9) 

#Response variables of interest; 120 total in master data
response_vars <- c("Bac_Shannon_soil","Copiotroph_rel_abund_soil", 
                   "Oligotroph_rel_abund_soil", "OtherBacteria_rel_abund_soil", 
                   "Fungal_Shannon_soil", "AMF_abund_soil", 
                   "ECM_abund_soil", "Endo_abund_soil", "Other_abund_soil", "Path_abund_soil", 
                   "Sap_abund_soil", "AMF_Shannon_soil",
                   "Bac_Shannon_root", "Copiotroph_rel_abund_root", 
                   "Oligotroph_rel_abund_root", "OtherBacteria_rel_abund_root", 
                   "Fungal_Shannon_root", "AMF_abund_root", 
                   "ECM_abund_root", "Endo_abund_root", "Other_abund_root", "Path_abund_root", 
                   "Sap_abund_root", "AMF_Shannon_root", "perc_N", 
                   "perc_C", "CN_ratio", "ergosterol", "tot_P", "perc_P", "soil_moisture", 
                   "organic_matter_LOI", "litter_depth", 
                   "pine_litter_prop", "euc_litter_prop",
                    "total_litter_biomass", "Pyruvic_Acid_Methyl_Ester_BiologDay5", 
                   "Tween_40_BiologDay5", "Tween_80_BiologDay5", "a.Cyclodextrin_BiologDay5", 
                   "Glycogen_BiologDay5", "D.Cellobiose_BiologDay5", "a.D.Lactose_BiologDay5", 
                   "B.Methyl.D.Glucoside_BiologDay5", "D.Xylose_BiologDay5", 
                   "i.Erythritol_BiologDay5", "D.Mannitol_BiologDay5", 
                   "N.Acetyl.D.Glucosamine_BiologDay5", "D.Glucosaminic_Acid_BiologDay5", 
                   "Glucose.1.Phosphae_BiologDay5", "D.L.a.Glycerol_Phosphate_BiologDay5", 
                   "D.Galactonic_Acid_gamma_Lactone_BiologDay5", "D.Galactonic_Acid_BiologDay5", 
                   "X2.Hydroxy_Benzoic_Acid_BiologDay5", "X4.Hydroxy_Benzoic_Acid_BiologDay5", 
                   "gamma.Amino_Butyric_Acid_BiologDay5", "Itaconic_Acid_BiologDay5", 
                   "a.Keto_Butyric_Acid_BiologDay5", "D.Malic_Acid_BiologDay5", 
                   "L.Arginine_BiologDay5", "L.Asparagine_BiologDay5", 
                   "L.Phenylalanine_BiologDay5", "L.Serine_BiologDay5", 
                   "L.Threonine_BiologDay5", "Glycyl.L.Glutamic_Acid_BiologDay5", 
                   "Phenylethyl.amine_BiologDay5", "Putrescine_BiologDay5", 
                   "Euc_leafLitter_ergosterol", "Euc_leafLitter_totP", "Euc_leafLitter_percP", 
                   "Euc_leafLitter_percN", "Euc_leafLitter_percC", "Pine_leafLitter_ergosterol", 
                   "Pine_leafLitter_totP", "Pine_leafLitter_percP", "Pine_leafLitter_percN", 
                   "Pine_leafLitter_percC", "Euc_Bark_ergosterol", "Euc_Bark_totP", 
                   "Euc_Bark_percP", "Euc_Bark_percN", "Euc_Bark_percC", 
                   "Litter_OLayer_ergosterol", "Litter_OLayer_totP", "Litter_OLayer_percP", 
                   "Litter_OLayer_percN", "Litter_OLayer_percC", "Litter_avg_ergosterol", 
                   "Litter_avg_totP", "Litter_avg_percP", "Litter_avg_percN", 
                   "Litter_avg_percC")

##########All response variables#################
# all_response_vars <- c("treeName", "tree_num", "tree_spp", "Bac_Shannon_soil", "Bac_Simpson_soil", 
#                    "Bac_Invsimpson_soil", "Bac_Richness_soil", "TreeSpecies", 
#                    "Copiotroph_rel_abund_soil", "Oligotroph_rel_abund_soil", 
#                    "OtherBacteria_rel_abund_soil", "Copiotroph_richness_soil", 
#                    "Oligotroph_richness_soil", "OtherBacteria_richness_soil", 
#                    "Fungal_Shannon_soil", "Fungal_Simpson_soil", "Fungal_Invsimpson_soil", 
#                    "Fungal_Richness_soil", "AMF_abund_soil", "ECM_abund_soil", "Endo_abund_soil", 
#                    "Other_abund_soil", "Path_abund_soil", "Sap_abund_soil", "AMF_richness_soil", 
#                    "ECM_richness_soil", "Endo_richness_soil", "Other_richness_soil", 
#                    "Path_richness_soil", "Sap_richness_soil", "AMF_Shannon_soil", 
#                    "AMF_Simpson_soil", "AMF_Invsimpson_soil", "AMF_Richness_soil", 
#                    "Bac_Shannon_root", "Bac_Simpson_root", "Bac_Invsimpson_root", 
#                    "Bac_Richness_root", "Copiotroph_rel_abund_root", "Oligotroph_rel_abund_root", 
#                    "OtherBacteria_rel_abund_root", "Copiotroph_richness_root", 
#                    "Oligotroph_richness_root", "OtherBacteria_richness_root", 
#                    "Fungal_Shannon_root", "Fungal_Simpson_root", "Fungal_Invsimpson_root", 
#                    "Fungal_Richness_root", "AMF_abund_root", "ECM_abund_root", "Endo_abund_root", 
#                    "Other_abund_root", "Path_abund_root", "Sap_abund_root", "AMF_richness_root", 
#                    "ECM_richness_root", "Endo_richness_root", "Other_richness_root", 
#                    "Path_richness_root", "Sap_richness_root", "AMF_Shannon_root", 
#                    "AMF_Simpson_root", "AMF_Invsimpson_root", "AMF_Richness_root", 
#                    "perc_N", "perc_C", "CN_ratio", "ergosterol", "tot_P", "perc_P", 
#                    "soil_moisture", "organic_matter_LOI", "dbh_cm", "height_m", 
#                    "litter_depth", "pine_litter_prop", "euc_litter_prop", "latitude", 
#                    "longitude", "lomandra_n", "lomandra_p", "lomandra_k", "lomandra_c", 
#                    "total_litter_biomass", "Water_BiologDay1", "Pyruvic_Acid_Methyl_Ester_BiologDay1", 
#                    "Tween_40_BiologDay1", "Tween_80_BiologDay1", "a.Cyclodextrin_BiologDay1", 
#                    "Glycogen_BiologDay1", "D.Cellobiose_BiologDay1", "a.D.Lactose_BiologDay1", 
#                    "B.Methyl.D.Glucoside_BiologDay1", "D.Xylose_BiologDay1", "i.Erythritol_BiologDay1", 
#                    "D.Mannitol_BiologDay1", "N.Acetyl.D.Glucosamine_BiologDay1", 
#                    "D.Glucosaminic_Acid_BiologDay1", "Glucose.1.Phosphae_BiologDay1", 
#                    "D.L.a.Glycerol_Phosphate_BiologDay1", "D.Galactonic_Acid_gamma_Lactone_BiologDay1", 
#                    "D.Galactonic_Acid_BiologDay1", "X2.Hydroxy_Benzoic_Acid_BiologDay1", 
#                    "X4.Hydroxy_Benzoic_Acid_BiologDay1", "gamma.Amino_Butyric_Acid_BiologDay1", 
#                    "Itaconic_Acid_BiologDay1", "a.Keto_Butyric_Acid_BiologDay1", 
#                    "D.Malic_Acid_BiologDay1", "L.Arginine_BiologDay1", "L.Asparagine_BiologDay1", 
#                    "L.Phenylalanine_BiologDay1", "L.Serine_BiologDay1", "L.Threonine_BiologDay1", 
#                    "Glycyl.L.Glutamic_Acid_BiologDay1", "Phenylethyl.amine_BiologDay1", 
#                    "Putrescine_BiologDay1", "Water_BiologDay3", "Pyruvic_Acid_Methyl_Ester_BiologDay3", 
#                    "Tween_40_BiologDay3", "Tween_80_BiologDay3", "a.Cyclodextrin_BiologDay3", 
#                    "Glycogen_BiologDay3", "D.Cellobiose_BiologDay3", "a.D.Lactose_BiologDay3", 
#                    "B.Methyl.D.Glucoside_BiologDay3", "D.Xylose_BiologDay3", "i.Erythritol_BiologDay3", 
#                    "D.Mannitol_BiologDay3", "N.Acetyl.D.Glucosamine_BiologDay3", 
#                    "D.Glucosaminic_Acid_BiologDay3", "Glucose.1.Phosphae_BiologDay3", 
#                    "D.L.a.Glycerol_Phosphate_BiologDay3", "D.Galactonic_Acid_gamma_Lactone_BiologDay3", 
#                    "D.Galactonic_Acid_BiologDay3", "X2.Hydroxy_Benzoic_Acid_BiologDay3", 
#                    "X4.Hydroxy_Benzoic_Acid_BiologDay3", "gamma.Amino_Butyric_Acid_BiologDay3", 
#                    "Itaconic_Acid_BiologDay3", "a.Keto_Butyric_Acid_BiologDay3", 
#                    "D.Malic_Acid_BiologDay3", "L.Arginine_BiologDay3", "L.Asparagine_BiologDay3", 
#                    "L.Phenylalanine_BiologDay3", "L.Serine_BiologDay3", "L.Threonine_BiologDay3", 
#                    "Glycyl.L.Glutamic_Acid_BiologDay3", "Phenylethyl.amine_BiologDay3", 
#                    "Putrescine_BiologDay3", "Water_BiologDay5", "Pyruvic_Acid_Methyl_Ester_BiologDay5", 
#                    "Tween_40_BiologDay5", "Tween_80_BiologDay5", "a.Cyclodextrin_BiologDay5", 
#                    "Glycogen_BiologDay5", "D.Cellobiose_BiologDay5", "a.D.Lactose_BiologDay5", 
#                    "B.Methyl.D.Glucoside_BiologDay5", "D.Xylose_BiologDay5", "i.Erythritol_BiologDay5", 
#                    "D.Mannitol_BiologDay5", "N.Acetyl.D.Glucosamine_BiologDay5", 
#                    "D.Glucosaminic_Acid_BiologDay5", "Glucose.1.Phosphae_BiologDay5", 
#                    "D.L.a.Glycerol_Phosphate_BiologDay5", "D.Galactonic_Acid_gamma_Lactone_BiologDay5", 
#                    "D.Galactonic_Acid_BiologDay5", "X2.Hydroxy_Benzoic_Acid_BiologDay5", 
#                    "X4.Hydroxy_Benzoic_Acid_BiologDay5", "gamma.Amino_Butyric_Acid_BiologDay5", 
#                    "Itaconic_Acid_BiologDay5", "a.Keto_Butyric_Acid_BiologDay5", 
#                    "D.Malic_Acid_BiologDay5", "L.Arginine_BiologDay5", "L.Asparagine_BiologDay5", 
#                    "L.Phenylalanine_BiologDay5", "L.Serine_BiologDay5", "L.Threonine_BiologDay5", 
#                    "Glycyl.L.Glutamic_Acid_BiologDay5", "Phenylethyl.amine_BiologDay5", 
#                    "Putrescine_BiologDay5", "Euc_leafLitter_ergosterol", "Euc_leafLitter_totP", 
#                    "Euc_leafLitter_per
# #############
#Initialize an empty list to store effect sizes
effectsize_results <- list()

#Loop through each response variable to calculate effect sizes using Cohen's d
for (var in response_vars) {
  #Calculate means, SDs, and sample sizes for both Eucalyptus & Pine
  mean_euc <- mean(AU_masterdata[[var]][AU_masterdata$tree_spp == "euc"], na.rm = TRUE)
  sd_euc <- sd(AU_masterdata[[var]][AU_masterdata$tree_spp == "euc"], na.rm = TRUE)
  n_euc <- sum(AU_masterdata$tree_spp == "euc", na.rm = TRUE)
  
  mean_pine <- mean(AU_masterdata[[var]][AU_masterdata$tree_spp == "pine"], na.rm = TRUE)
  sd_pine <- sd(AU_masterdata[[var]][AU_masterdata$tree_spp == "pine"], na.rm = TRUE)
  n_pine <- sum(AU_masterdata$tree_spp == "pine", na.rm = TRUE)
  
  #Calculate Cohen's d effect size for Pine vs Eucalyptus
  effect_size <- escalc(measure = "SMD", m1i = mean_pine, m2i = mean_euc,
                        sd1i = sd_pine, sd2i = sd_euc,
                        n1i = n_pine, n2i = n_euc)
  
  #Store results in the list w/ variable name as row name
  effectsize_results[[var]] <- data.frame(
    Variable = var,
    EffectSize = effect_size$yi,
    Variance = effect_size$vi,
    AbsEffectSize = abs(effect_size$yi) # for sorting by largest effect size
  )
}

#Combine results into a single data frame
effect_data <- bind_rows(effectsize_results)

#Select the top [ex: 10] response variables with the largest absolute effect sizes
top_effects <- effect_data %>%
  arrange(desc(AbsEffectSize)) %>%
  slice(1:10) #change depending on how many top effect size values you want to display

#Meta-analysis model and forest plot for the top [ex: 10] variables
meta_model <- rma(yi = EffectSize, vi = Variance, data = top_effects)

#Plot the forest plot of effect sizes for the top 10 response variables
forest(meta_model, slab = top_effects$Variable, header = c("Response Variable", "Effect Size (Pine vs. Euc)"))
#####################


########Cohen's d, response variables grouped######
AU_masterdata <- read.csv("/Users/elenaleander/Desktop/AU_manuscript/Aus23_allData_19Nov24.csv", header = TRUE)

#Define response variable groups
ecoplate_vars <- c("Pyruvic_Acid_Methyl_Ester_BiologDay5", 
                   "Tween_40_BiologDay5", "Tween_80_BiologDay5", "a.Cyclodextrin_BiologDay5", 
                   "Glycogen_BiologDay5", "D.Cellobiose_BiologDay5", "a.D.Lactose_BiologDay5", 
                   "B.Methyl.D.Glucoside_BiologDay5", "D.Xylose_BiologDay5", 
                   "i.Erythritol_BiologDay5", "D.Mannitol_BiologDay5", 
                   "N.Acetyl.D.Glucosamine_BiologDay5", "D.Glucosaminic_Acid_BiologDay5", 
                   "Glucose.1.Phosphae_BiologDay5", "D.L.a.Glycerol_Phosphate_BiologDay5", 
                   "D.Galactonic_Acid_gamma_Lactone_BiologDay5", "D.Galactonic_Acid_BiologDay5", 
                   "X2.Hydroxy_Benzoic_Acid_BiologDay5", "X4.Hydroxy_Benzoic_Acid_BiologDay5", 
                   "gamma.Amino_Butyric_Acid_BiologDay5", "Itaconic_Acid_BiologDay5", 
                   "a.Keto_Butyric_Acid_BiologDay5", "D.Malic_Acid_BiologDay5", 
                   "L.Arginine_BiologDay5", "L.Asparagine_BiologDay5", 
                   "L.Phenylalanine_BiologDay5", "L.Serine_BiologDay5", 
                   "L.Threonine_BiologDay5", "Glycyl.L.Glutamic_Acid_BiologDay5", 
                   "Phenylethyl.amine_BiologDay5", "Putrescine_BiologDay5")

soil_microbial_vars <- c("Bac_Shannon_soil","Copiotroph_rel_abund_soil", 
                         "Oligotroph_rel_abund_soil", "OtherBacteria_rel_abund_soil", 
                         "Fungal_Shannon_soil", "AMF_abund_soil", 
                         "ECM_abund_soil", "Endo_abund_soil", "Other_abund_soil", "Path_abund_soil", 
                         "Sap_abund_soil", "AMF_Shannon_soil")

root_microbial_vars <- c("Bac_Shannon_root", "Copiotroph_rel_abund_root", 
                        "Oligotroph_rel_abund_root", "OtherBacteria_rel_abund_root", 
                        "Fungal_Shannon_root", "AMF_abund_root", 
                        "ECM_abund_root", "Endo_abund_root", "Other_abund_root", "Path_abund_root", 
                        "Sap_abund_root", "AMF_Shannon_root")

litter_vars <- c("Euc_leafLitter_ergosterol", "Euc_leafLitter_totP", "Euc_leafLitter_percP", 
                 "Euc_leafLitter_percN", "Euc_leafLitter_percC", "Pine_leafLitter_ergosterol", 
                 "Pine_leafLitter_totP", "Pine_leafLitter_percP", "Pine_leafLitter_percN", 
                 "Pine_leafLitter_percC", "Euc_Bark_ergosterol", "Euc_Bark_totP", 
                 "Euc_Bark_percP", "Euc_Bark_percN", "Euc_Bark_percC", 
                 "Litter_OLayer_ergosterol", "Litter_OLayer_totP", "Litter_OLayer_percP", 
                 "Litter_OLayer_percN", "Litter_OLayer_percC", "Litter_avg_ergosterol", 
                 "Litter_avg_totP", "Litter_avg_percP", "Litter_avg_percN", 
                 "Litter_avg_percC","litter_depth", 
                 "pine_litter_prop", "euc_litter_prop",
                 "total_litter_biomass") 

soil_vars <- c("perc_N","perc_C", "CN_ratio", "ergosterol", "tot_P", "perc_P", "soil_moisture", 
                 "organic_matter_LOI")

#Calculate effect sizes for a group of variables and return a data frame
calculate_effect_sizes <- function(variables, data) {
  effectsize_results <- list()
  
  for (var in variables) {
    #Calculate means, SDs, and sample sizes for Eucalyptus & Pine
    mean_euc <- mean(data[[var]][data$tree_spp == "euc"], na.rm = TRUE)
    sd_euc <- sd(data[[var]][data$tree_spp == "euc"], na.rm = TRUE)
    n_euc <- sum(data$tree_spp == "euc", na.rm = TRUE)
    
    mean_pine <- mean(data[[var]][data$tree_spp == "pine"], na.rm = TRUE)
    sd_pine <- sd(data[[var]][data$tree_spp == "pine"], na.rm = TRUE)
    n_pine <- sum(data$tree_spp == "pine", na.rm = TRUE)
    
    #Calculate Cohen's d effect size for Pine vs Eucalyptus
    effect_size <- escalc(measure = "SMD", m1i = mean_pine, m2i = mean_euc,
                          sd1i = sd_pine, sd2i = sd_euc,
                          n1i = n_pine, n2i = n_euc)
    
    #Calculate confidence intervals to determine significance
    ci_low <- effect_size$yi - 1.96 * sqrt(effect_size$vi)
    ci_high <- effect_size$yi + 1.96 * sqrt(effect_size$vi)
    significant <- ifelse(ci_low > 0 | ci_high < 0, "Significant", "Not Significant")
    
    #Store results in a data frame
    effectsize_results[[var]] <- data.frame(
      Variable = var,
      EffectSize = effect_size$yi,
      Variance = effect_size$vi,
      CI_Lower = ci_low,
      CI_Upper = ci_high,
      Significant = significant  # Add significance column
    )
  }
  
  #Combine results into a single data frame
  effect_data <- bind_rows(effectsize_results)
  return(effect_data)
}

#Calculate effect sizes for each group with significance flags
ecoplate_effects <- calculate_effect_sizes(ecoplate_vars, AU_masterdata)
soil_microbial_effects <- calculate_effect_sizes(soil_microbial_vars, AU_masterdata)
root_microbial_effects <- calculate_effect_sizes(root_microbial_vars, AU_masterdata)
litter_effects <- calculate_effect_sizes(litter_vars, AU_masterdata)
soil_effects <- calculate_effect_sizes(soil_vars, AU_masterdata)

#Create a forest plot for a given set of effect sizes
create_forest_plot <- function(effect_data, title) {
  #Run the meta-analysis model
  meta_model <- rma(yi = EffectSize, vi = Variance, data = effect_data)
  
  #Create the forest plot
  forest(meta_model, slab = effect_data$Variable, header = c("Response Variable", "Effect Size (Pine vs. Euc)"))
  
  #Add title
  title(main = title)
}

#Generate forest plots for each group
create_forest_plot(ecoplate_effects, "EcoPlate Response Variables")
create_forest_plot(soil_microbial_effects, "Soil Microbial Response Variables")
create_forest_plot(root_microbial_effects, "Root Microbial Response Variables")
create_forest_plot(litter_effects, "Litter Properties Response Variables")
create_forest_plot(soil_effects, "Soil Properties Response Variables")

#Create a table of significant effect sizes for a given group
create_significant_effects_table <- function(effect_data, group_name) {
  # Filter only significant effects
  significant_effects <- effect_data %>%
    filter(Significant == "Significant") %>%
    mutate(Group = group_name)  # Add a column for the group name
  return(significant_effects)
}

#Apply this to each group of response variables
ecoplate_significant <- create_significant_effects_table(ecoplate_effects, "EcoPlate")
soil_microbial_significant <- create_significant_effects_table(soil_microbial_effects, "Soil Microbial")
root_microbial_significant <- create_significant_effects_table(root_microbial_effects, "Root Microbial")
litter_significant <- create_significant_effects_table(litter_effects, "Litter")
soil_significant <- create_significant_effects_table(soil_effects, "Soil")

#Combine all significant effect tables into one data frame
all_significant_effects <- bind_rows(
  ecoplate_significant, soil_microbial_significant, 
  root_microbial_significant, litter_significant, soil_significant
)

#Select and arrange columns for the final table
significant_effects_table <- all_significant_effects %>%
  select(Group, Variable, EffectSize, CI_Lower, CI_Upper) %>%  # Choose desired columns
  arrange(Group, desc(EffectSize))  # Sort by group and effect size

#Display the table of significant response variables and their effect sizes
print(significant_effects_table)

##########Cohen's d, only significant response variables#########
#Combine all significant effects into a single data frame [*see above]
all_significant_effects <- bind_rows(
  ecoplate_significant, soil_microbial_significant, 
  root_microbial_significant, litter_significant, soil_significant
)

#Filter the data to only significant variables
significant_only <- all_significant_effects %>%
  select(Group, Variable, EffectSize, CI_Lower, CI_Upper, Variance) %>%
  arrange(Group, desc(EffectSize))

#Create forest plot of significant response variables only
create_significant_forest_plot <- function(effect_data) {
  meta_model <- rma(yi = EffectSize, vi = Variance, data = effect_data)
  forest(meta_model, slab = paste(effect_data$Group, effect_data$Variable, sep = " - "),
         header = c("Response Variable", "Effect Size (Pine vs. Euc)"))
  title(main = "Significant Response Variables")
}

create_significant_forest_plot(significant_only)

#Add clean labels
create_significant_forest_plot <- function(effect_data) {

  custom_labels <- c(
    "EcoPlate - Tween 40",
    "EcoPlate - L-serine",
    "EcoPlate - Glycyl-L-glutamic acid",
    "EcoPlate - Erythritol",
    "EcoPlate - Putrescine",
    "EcoPlate - L-putrescine",
    "Litter - Proportion Pine",
    "Litter - Eucalyptus % P",
    "Litter - Total P Eucalyptus",
    "Litter - O Layer Ergosterol",
    "Litter - Average Ergosterol",
    "Litter - Proportion Eucalyptus",
    "Root Microbial - Pathogenic Abundance",
    "Soil - Ergosterol",
    "Soil - Soil Moisture",
    "Soil Microbial - Bacteria Shannon Diversity"
  )
  
  #Fit meta-analytic model
  meta_model <- rma(yi = EffectSize, vi = Variance, data = effect_data)
  
  #Create forest plot with custom labels
  forest(meta_model, slab = custom_labels, 
         header = c("Response Variable", "Effect Size (Pine vs. Euc)"))
  

  title(main = "Significant Response Variables")
}

create_significant_forest_plot(significant_only)

