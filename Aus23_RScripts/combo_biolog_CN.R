#BB
#10/23/24
#Merging the pooled soil data from Jason with the ecoplate data generated from Elena/Edith.


library(dplyr)
library(tidyr)


eco <- read.csv("Merged_data/eco_plates_wide_averaged.csv")
cn <- read.csv("Merged_data/Aus23_microbes_soilCNP_pooled.csv")

cn$treeName

eco <- eco %>%
  mutate(treeName = case_when(
    tree_spp == "euc" ~ paste0("E", tree_num),
    tree_spp == "pine" ~ paste0("P", tree_num),
    TRUE ~ NA_character_  # Handles other cases if needed
  ))

ds <- full_join(eco, cn)

write.csv(ds, "Merged_data/Aus23_CNP_pooled_biolog.csv")
