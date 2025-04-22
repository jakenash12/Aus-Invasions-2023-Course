# Look at the most abundant OTUs in each sample for pine soil
# 4/17/25 CV

library(tidyverse)

# read in rarefied OTU tables
otu.its =
  read_delim("Aus-Invasions-2023-Course/Aus23_ITS_Metabarcoding/OTUTables/OtuMatITS_rare_pooled.csv", delim=",") %>%
  column_to_rownames(var = colnames(.)[1])
# read in taxonomy
its.tax <- read_tsv("Aus-Invasions-2023-Course/Aus23_ITS_Metabarcoding/ITS2_Dada2_repseqs97_taxonomy_edited_FunT_FunG.tsv")
its.tax.split <- its.tax %>% 
  separate(Kingdom, into = c("Phylum", "Class", "Order", "Family", "Genus"), sep = ";") %>%
  column_to_rownames(var = "Feature ID") %>%
  mutate(Kingdom = "Fungi") %>%
  filter(str_detect(guild, "Ectomycorrhizal")) # keep only EMF taxa
  select(Species) %>%
  mutate(across(everything(), ~na_if(., "NA"))) %>%
  mutate(across(everything(), ~ifelse(grepl("Incertae_sedis", ., fixed = TRUE), NA, .)))

# adjust so OTU names are rows
otu.its.t <- as.data.frame(t(otu.its))
otu.its.t <- merge(otu.its.t, its.tax.split, by = 0, all.x=TRUE)
colnames(otu.its.t)[1] <- "OTU"
otu.its.t <- otu.its.t[! is.na(otu.its.t$Species),]

# select only pine soil samples
pine_cols <- grep("P", names(otu.its.t), value = TRUE)
pine_soil_cols <- pine_cols[grepl("Soil", pine_cols)]
pine_soil_cols2 <- c(pine_soil_cols, "Species")
pine_soil <- otu.its.t[, pine_soil_cols2]

# Sort to see most abundant EMF taxa in each pine soil sample

### >40cm DBH ###
# 70cm DBH
P1 <- pine_soil[c("Species", "P1_Soil")]
P1_sorted <- P1 %>% arrange(desc(P1_Soil))

# 84cm DBH
P5 <- pine_soil[c("Species", "P5_Soil")]
P5_sorted <- P5 %>% arrange(desc(P5_Soil))

# 43cm DBH
P8 <- pine_soil[c("Species", "P8_Soil")]
P8_sorted <- P8 %>% arrange(desc(P8_Soil))

# 45cm DBH
P9 <- pine_soil[c("Species", "P9_Soil")]
P9_sorted <- P9 %>% arrange(desc(P9_Soil))


### <40cm DBH ###

# 16cm DBH
P2 <- pine_soil[c("Species", "P2_Soil")]
P2_sorted <- P2 %>% arrange(desc(P2_Soil))

# 27cm DBH
P3 <- pine_soil[c("Species", "P3_Soil")]
P3_sorted <- P3 %>% arrange(desc(P3_Soil))

# 29cm DBH
P4 <- pine_soil[c("Species", "P4_Soil")]
P4_sorted <- P4 %>% arrange(desc(P4_Soil))

# 37cm DBH
P6 <- pine_soil[c("Species", "P6_Soil")]
P6_sorted <- P6 %>% arrange(desc(P6_Soil))

# 22cm DBH
P7 <- pine_soil[c("Species", "P7_Soil")]
P7_sorted <- P7 %>% arrange(desc(P7_Soil))





