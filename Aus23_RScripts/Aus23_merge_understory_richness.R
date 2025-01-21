# Merging master dataframe with understory species richness data
# 1/18/24 C Vietorisz

library(tidyverse)

setwd("Aus-Invasions-2023-Course/")

# read in main merged dataframe
main <- read.csv("Merged_data/Aus23_allData_19Nov24.csv", row.names=1)

# read in understory species richness data
rich <- read.csv("Plant_data/Belanglo_understory_richness.csv")

# sum the N and S richness for each tree to get plant morphotype richness per tree
rich.sum <- rich %>%
  group_by(Focal_tree_ID) %>%
  summarise(Understory_richness_perTree = sum(Sp_richness_adj)) %>%
  rename(treeName = Focal_tree_ID)

# merge all data together
main.rich <- merge(main, rich.sum, by = "treeName")

write.csv(main.rich, "Merged_data/Aus23_allData_18Jan25.csv")
