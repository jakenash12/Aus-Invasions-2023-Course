# Make a heatmap of differentially expressed OTUs under pine vs euc against environmental variables
# 11/21/24
# CV

library(DescTools)

setwd("~/Desktop/GitHub/Aus-Invasions-2023-Course")

# read in 16S OTU table (converted to relative abundance)
otu16s <- read.csv("Aus23_16S_Metabarcoding/OTUTables/OtuMat16S_rel.csv", row.names=1)
otu16s$treeName <- substr(rownames(otu16s), 1, 2)
otu16s$SampleType <- sub(".*_", "", rownames(otu16s))

# subset out soil and roots separately
otu16s.root <- subset(otu16s, otu16s$SampleType %like% "Root")
otu16s.soil <- subset(otu16s, otu16s$SampleType %like% "Soil")

# read in differentially expressed OTUs
root.de.otu <- 

# read in full data
fulldat <- read.csv("Merged_data/Aus23_allData_19Nov24.csv", row.names=1)
sigvars <- c("Tween_40_BiologDay5", "L.Serine_BiologDay5", "Glycyl.L.Glutamic_Acid_BiologDay5",
               "i.Erythritol_BiologDay5", "Putrescine_BiologDay5", "L.Asparagine_BiologDay5", 
               "pine_litter_prop", "Euc_leafLitter_percP", "Euc_leafLitter_totP",              
               "Litter_OLayer_ergosterol", "Litter_avg_ergosterol","euc_litter_prop", "Path_abund_root", "ergosterol","soil_moisture",
               "Bac_Shannon_soil")


