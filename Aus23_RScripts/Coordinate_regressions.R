# Regressions incorporating latitude and longitude on the abundances of genera enriched in pines or eucalypts
# 4/22/25 CRV

library(car)
library(emmeans)
library(tidyverse)

setwd("Aus-Invasions-2023-Course/")

# read in genus-level pooled OTU tables
bac.root.gen = read_delim("Aus23_16S_Metabarcoding/OTUTables/OtuMat16S_rare_Roots_Top5Genus_pooled.csv") %>%
  column_to_rownames(var = colnames(.)[1])
bac.soil.gen = read_delim("Aus23_16S_Metabarcoding/OTUTables/OtuMat16S_rare_Soil_Top5Genus_pooled.csv") %>%
  column_to_rownames(var = colnames(.)[1])
fun.root.gen = read_delim("Aus23_ITS_Metabarcoding/OTUTables/OtuMatITS_rare_Roots_Top5Genus_pooled.csv") %>%
  column_to_rownames(var = colnames(.)[1])
fun.soil.gen = read_delim("Aus23_ITS_Metabarcoding/OTUTables/OtuMatITS_rare_Soil_Top5Genus_pooled.csv") %>%
  column_to_rownames(var = colnames(.)[1])

# add with coordinates, tree species, and DBH
coords <- read.csv("Merged_data/AUlatlong.csv")
dbh <- read.csv("Merged_data/Aus23_allData_18Jan25.csv")[c("treeName", "dbh_cm")]
data <- merge(coords, dbh, by = "treeName")

# create dataframes by sample type
bac.root.data <- merge(bac.root.gen, data, by = "treeName")
bac.soil.data <- merge(bac.soil.gen, data, by = "treeName")
fun.root.data <- merge(fun.root.gen, data, by = "treeName")
fun.soil.data <- merge(fun.soil.gen, data, by = "treeName")

# move treeName to the end of the dataframe so the genera abundances are in columns 1:10
bac.root.data <- bac.root.data %>% relocate(treeName, .after = last_col())
bac.soil.data <- bac.soil.data %>% relocate(treeName, .after = last_col())
fun.root.data <- fun.root.data %>% relocate(treeName, .after = last_col())
fun.soil.data <- fun.soil.data %>% relocate(treeName, .after = last_col())

# FUNGAL ROOTS
# create empty matrix to hold the stats
fun.root.stats <- matrix(nrow=10, ncol=10)
colnames(fun.root.stats) <- c("Latitutde_F", "Longitude_F", "DBH_F", "tree_spp_F", "DBHxTree_F", "Latitutde_P", "Longitude_P", "DBH_P", "tree_spp_P", "DBHxTree_P")
rownames(fun.root.stats) <- colnames(fun.root.data)[1:10]
# run a for loop to run the lm on all genera
for (i in 1:10) {
lm.result=lm(fun.root.data[[i]] ~ latitude	+ longitude+dbh_cm*tree_spp, data = fun.root.data)
lm.anova <- Anova(lm.result)
F_vals <- lm.anova$`F value`[1:5]
P_vals <- lm.anova$`Pr(>F)`[1:5]
fun.root.stats[i,1:5] <- F_vals
fun.root.stats[i,6:10] <- P_vals
}

# FUNGAL SOILS
# create empty matrix to hold the stats
fun.soil.stats <- matrix(nrow=10, ncol=10)
colnames(fun.soil.stats) <- c("Latitutde_F", "Longitude_F", "DBH_F", "tree_spp_F", "DBHxTree_F", "Latitutde_P", "Longitude_P", "DBH_P", "tree_spp_P", "DBHxTree_P")
rownames(fun.soil.stats) <- colnames(fun.soil.data)[1:10]
# run a for loop to run the lm on all genera
for (i in 1:10) {
  lm.result=lm(fun.soil.data[[i]] ~ latitude	+ longitude+dbh_cm*tree_spp, data = fun.soil.data)
  lm.anova <- Anova(lm.result)
  F_vals <- lm.anova$`F value`[1:5]
  P_vals <- lm.anova$`Pr(>F)`[1:5]
  fun.soil.stats[i,1:5] <- F_vals
  fun.soil.stats[i,6:10] <- P_vals
}

# BACTERIAL ROOTS
# create empty matrix to hold the stats
bac.root.stats <- matrix(nrow=10, ncol=10)
colnames(bac.root.stats) <- c("Latitutde_F", "Longitude_F", "DBH_F", "tree_spp_F", "DBHxTree_F", "Latitutde_P", "Longitude_P", "DBH_P", "tree_spp_P", "DBHxTree_P")
rownames(bac.root.stats) <- colnames(bac.root.data)[1:10]
# run a for loop to run the lm on all genera
for (i in 1:10) {
  lm.result=lm(bac.root.data[[i]] ~ latitude	+ longitude+dbh_cm*tree_spp, data = bac.root.data)
  lm.anova <- Anova(lm.result)
  F_vals <- lm.anova$`F value`[1:5]
  P_vals <- lm.anova$`Pr(>F)`[1:5]
  bac.root.stats[i,1:5] <- F_vals
  bac.root.stats[i,6:10] <- P_vals
}

# BACTERIAL SOILS
# create empty matrix to hold the stats
bac.soil.stats <- matrix(nrow=10, ncol=10)
colnames(bac.soil.stats) <- c("Latitutde_F", "Longitude_F", "DBH_F", "tree_spp_F", "DBHxTree_F", "Latitutde_P", "Longitude_P", "DBH_P", "tree_spp_P", "DBHxTree_P")
rownames(bac.soil.stats) <- colnames(bac.soil.data)[1:10]
# run a for loop to run the lm on all genera
for (i in 1:10) {
  lm.result=lm(bac.soil.data[[i]] ~ latitude	+ longitude+dbh_cm*tree_spp, data = bac.soil.data)
  lm.anova <- Anova(lm.result)
  F_vals <- lm.anova$`F value`[1:5]
  P_vals <- lm.anova$`Pr(>F)`[1:5]
  bac.soil.stats[i,1:5] <- F_vals
  bac.soil.stats[i,6:10] <- P_vals
}

# save outputs
write.csv(bac.soil.stats, "Autocorrelation_analyses/bac_soil_coordinate_regressions.csv")
write.csv(bac.root.stats, "Autocorrelation_analyses/bac_root_coordinate_regressions.csv")
write.csv(fun.soil.stats, "Autocorrelation_analyses/fungal_soil_coordinate_regressions.csv")
write.csv(fun.root.stats, "Autocorrelation_analyses/fungal_root_coordinate_regressions.csv")


# Bacterial Shannon diversity
shan <- read.csv("Merged_data/Aus23_allData_18Jan25.csv")#[c("treeName", "dbh_cm", "Bac_Shannon_soil")]
shan_full <- merge(coords, data2, by = "treeName")
lm.result=lm(Bac_Shannon_soil ~ latitude	+ longitude+dbh_cm*tree_spp, data = shan_full)
Anova(lm.result)

