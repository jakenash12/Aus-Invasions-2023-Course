# Differentially abundant genus heatmap
# 1/6/24 CV

library(DescTools)
library(tidyverse)

# read in full data
fulldat <- read.csv("Aus-Invasions-2023-Course/Merged_data/Aus23_allData_19Nov24.csv", row.names=1)

# read in genus-level pooled OTU tables
otu.genus.16S <- read.csv("Aus-Invasions-2023-Course/Aus23_16S_Metabarcoding/OTUTables/OtuMat16S_rel_byGenus_pooled.csv", row.names=1)
otu.genus.its <- read.csv("Aus-Invasions-2023-Course/Aus23_ITS_Metabarcoding/OTUTables/OtuMatITS_rel_byGenus_pooled.csv", row.names=1)
# add tree ID
otu.genus.16S <- otu.genus.16S %>%
  rownames_to_column("SampleID") %>%
  mutate(treeName = str_sub(SampleID, 1, 2))
otu.genus.its <- otu.genus.its %>%
  rownames_to_column("SampleID") %>%
  mutate(treeName = str_sub(SampleID, 1, 2))
  
# separate out by root and soil
soil.16s <- subset(otu.genus.16S, otu.genus.16S$SampleID %like% "%Soil")
root.16s <- subset(otu.genus.16S, otu.genus.16S$SampleID %like% "%Root")
soil.its <- subset(otu.genus.its, otu.genus.its$SampleID %like% "%Soil")
root.its <- subset(otu.genus.its, otu.genus.its$SampleID %like% "%Root")

# subset out for the genera we want
bac.root.gen.names <- c("treeName", "Lechevalieria", "Nocardia", "Acidicaldus", "SBR1031", "Novosphingobium", "Collimonas", "A4b")
fun.soil.gen.names <- c("treeName", "Phialocephala", "Thelephora", "Athelopsis", "Tomentella", "Leohumicola", "Clavulinopsis")
# why is there no Thelephora or Athelopsis in the data???

bac.root.gen <- root.16s[bac.root.gen.names]
fun.soil.gen <- soil.its[fun.soil.gen.names]

# merge each OTU table with environmental data



# pull out significant environmental variables and differentially abundant genera
sigvars <- c("treeName", "Tween_40_BiologDay5", "L.Serine_BiologDay5", "Glycyl.L.Glutamic_Acid_BiologDay5",
             "i.Erythritol_BiologDay5", "Putrescine_BiologDay5", "L.Asparagine_BiologDay5", 
             "pine_litter_prop", "Euc_leafLitter_percP", "Euc_leafLitter_totP",              
             "Litter_OLayer_ergosterol", "Litter_avg_ergosterol","euc_litter_prop", "Path_abund_root", "ergosterol","soil_moisture",
             "Bac_Shannon_soil")
sigdat <- fulldat[sigvars]





