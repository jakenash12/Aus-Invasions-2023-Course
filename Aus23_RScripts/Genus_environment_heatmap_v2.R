# Differentially abundant genus heatmaps
# 1/18/24 CV
# removing p value stars from heatmaps

library(DescTools)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("Aus-Invasions-2023-Course/")

# read in full data
fulldat <- read.csv("Merged_data/Aus23_allData_18Jan25.csv", row.names=1)

# read in genus-level pooled OTU tables
otu.genus.16S = read_delim("Aus23_16S_Metabarcoding/OTUTables/OtuMat16S_rare_byGenus_pooled.csv") %>%
  column_to_rownames(var = colnames(.)[1])
otu.genus.its = read_delim("Aus23_ITS_Metabarcoding/OTUTables/OtuMatITS_rare_byGenus_pooled.csv") %>%
  column_to_rownames(var = colnames(.)[1])

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
bac.root.gen.names <- c("treeName", "Lechevalieria", "Nocardia", "Acidicaldus", "Nevskia", "Occallatibacter", "SBR1031", "Novosphingobium", "Collimonas", "A4b", "Bryocella")
bac.soil.gen.names <- c("treeName", "Gryllotalpicola", "Limnobacter", "Acidisoma", "Dinghuibacter", "Acidiphilium", "11-24", "SM1A02", "Subgroup_5", "A4b", "Mesorhizobium")
fun.root.gen.names <- c("treeName", "Phialocephala", "Talaromyces", "Devriesia", "Aspergillus", "Penicillium", "Auricularia", "Pezoloma", "Cryptosporiopsis", "Pezicula", "Mollisia")
fun.soil.gen.names <- c("treeName", "Phialocephala", "Thelephora", "Athelopsis", "Diplodia", "Ascocorticium", "Tomentella", "Leohumicola", "Clavulinopsis", "Pseudocamarosporium", "Flavocillium")

bac.root.gen <- root.16s[bac.root.gen.names]
bac.soil.gen <- soil.16s[bac.soil.gen.names]
fun.root.gen <- root.its[fun.root.gen.names]
fun.soil.gen <- soil.its[fun.soil.gen.names]

# pull out significant environmental variables and differentially abundant genera
sigvars <- c("treeName", "Tween_40_BiologDay5", "L.Serine_BiologDay5", "Glycyl.L.Glutamic_Acid_BiologDay5",
             "i.Erythritol_BiologDay5", "Putrescine_BiologDay5", "L.Asparagine_BiologDay5", 
             "pine_litter_prop", "Euc_leafLitter_percP",           
             "Litter_OLayer_ergosterol", "Litter_avg_ergosterol","euc_litter_prop", "ergosterol","soil_moisture",
             "Bac_Shannon_soil")
sigvar.order <-c("treeName", "Tween_40_BiologDay5", "L.Serine_BiologDay5", 
    "Glycyl.L.Glutamic_Acid_BiologDay5",
    "i.Erythritol_BiologDay5", "Putrescine_BiologDay5", "L.Asparagine_BiologDay5", "soil_moisture", "ergosterol", 
    "Litter_OLayer_ergosterol", "Litter_avg_ergosterol",
    "pine_litter_prop", "euc_litter_prop", "Euc_leafLitter_percP", "Bac_Shannon_soil")
sigdat <- fulldat[sigvars]
sigdat <- sigdat[,sigvar.order]
sigvar.names <- c("Tween 40 utilization", "L-Serine utilization", "Glycyl-L-Glutamic Acid utilization",
  "i-Erythritol utilization", "Putrescine utilization", "L-Asparagine utilization", "Soil moisture",
  "Soil ergosterol", "Litter O Layer ergosterol", "All litter layer ergosterol",
  "Proportion pine litter", "Proportion eucalypt litter", "Eucalypt leaf litter %P", "Bacterial Shannon diversity")
colnames(sigdat)[2:15] <- sigvar.names
                            
# merge each OTU table with environmental data
bac.root.sig <- merge(bac.root.gen, sigdat, by = "treeName")
bac.soil.sig <- merge(bac.soil.gen, sigdat, by = "treeName")
fun.root.sig <- merge(fun.root.gen, sigdat, by = "treeName")
fun.soil.sig <- merge(fun.soil.gen, sigdat, by = "treeName")

##################### Run correlations ##################### 
library(Hmisc)

X.order <-c("Tween 40 utilization", "L-Serine utilization", "Glycyl-L-Glutamic Acid utilization",
            "i-Erythritol utilization", "Putrescine utilization", "L-Asparagine utilization", "Soil moisture",
            "Soil ergosterol", "Litter O Layer ergosterol", "All litter layer ergosterol",
            "Proportion pine litter", "Proportion eucalypt litter", "Eucalypt leaf litter %P")

#### Bacterial roots #####

# Get correlation value (R) 
bac.root.r <- rcorr(as.matrix(bac.root.sig[2:ncol(bac.root.sig)]))$r
bac.root.r2 <- rcorr(as.matrix(bac.root.sig[2:ncol(bac.root.sig)]), type = "spearman")$r

# subset for just the relationships we want
bac.root.r.format <- bac.root.r %>% 
  as.data.frame(bac.root.r)  %>%
  dplyr::select(c(1:(length(bac.root.gen.names)-1))) %>% 
  slice(-(c(1:(length(bac.root.gen.names)-1), length(bac.root.r[,1]))))
rownames(bac.root.r.format) <- sigvar.names[1:13]

# pivot df longer
bac.root_long <- bac.root.r.format %>%
  rownames_to_column("Variable") %>%
  pivot_longer(cols = -Variable, names_to = "Genus", values_to = "Corr_value")

# make heatmap
bac.root.Yorder <- c("Lechevalieria", "Nocardia", "Acidicaldus", "Nevskia", "Occallatibacter", "SBR1031", "Novosphingobium", "Collimonas", "A4b", "Bryocella")
ggplot(bac.root_long, aes(x = factor(Variable, X.order), y = factor(Genus, bac.root.Yorder), fill = Corr_value)) +
  geom_tile() +
  scale_fill_gradient2(low = "dodgerblue4", mid = "white", high = "green4", midpoint = 0) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(fill = "Correlation") +
  coord_fixed()
ggsave("Correlations/Heatmaps/bacteria_root_genus_heatmap_v2.png", width = 5.5, heigh = 3.5, dpi = 300)

#### Bacterial soil #####

# Get correlation value (R) and P value (P)
bac.soil.r <- rcorr(as.matrix(bac.soil.sig[2:ncol(bac.soil.sig)]))$r

# subset for just the relationships we want
bac.soil.r.format <- bac.soil.r %>% 
  as.data.frame(bac.soil.r)  %>%
  dplyr::select(c(1:(length(bac.soil.gen.names)-1), length(bac.soil.r[,1]))) %>% 
  slice(-(c(1:(length(bac.soil.gen.names)-1), length(bac.soil.r[,1]))))
rownames(bac.soil.r.format) <- sigvar.names[1:13]

# pivot df longer
bac.soil_long <- bac.soil.r.format %>%
  rownames_to_column("Variable") %>%
  pivot_longer(cols = -Variable, names_to = "Genus", values_to = "Corr_value")

# make heatmap
bac.soil.Yorder <- c("Gryllotalpicola", "Limnobacter", "Acidisoma", "Dinghuibacter", "Acidiphilium", "11-24", "SM1A02", "Subgroup_5", "A4b", "Mesorhizobium", "Bacterial Shannon diversity")
ggplot(bac.soil_long, aes(x = factor(Variable, X.order), y = factor(Genus, bac.soil.Yorder), fill = Corr_value)) +
  geom_tile() +
  scale_fill_gradient2(low = "dodgerblue4", mid = "white", high = "green4", midpoint = 0) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(fill = "Correlation") +
  coord_fixed()
ggsave("Correlations/Heatmaps/bacteria_soil_genus_heatmap_v2.png", width = 5.5, heigh = 3.5, dpi = 300)


#### Fungal roots #####

# Get correlation value (R) and P value (P)
fun.root.r <- rcorr(as.matrix(fun.root.sig[2:ncol(fun.root.sig)]))$r

# subset for just the relationships we want
fun.root.r.format <- fun.root.r %>% 
  as.data.frame(fun.root.r)  %>%
  dplyr::select(1:(length(fun.root.gen.names)-1)) %>%
  slice(-(c(1:(length(fun.root.gen.names)-1), length(fun.root.r[,1]))))
rownames(fun.root.r.format) <- sigvar.names[1:13]

# pivot df longer
fun.root_long <- fun.root.r.format %>%
  rownames_to_column("Variable") %>%
  pivot_longer(cols = -Variable, names_to = "Genus", values_to = "Corr_value")

# make heatmap
fun.root.Yorder <- c("Phialocephala", "Talaromyces", "Devriesia", "Aspergillus", "Penicillium", "Auricularia", "Pezoloma", "Cryptosporiopsis", "Pezicula", "Mollisia")
ggplot(fun.root_long, aes(x = factor(Variable, X.order), y = factor(Genus, fun.root.Yorder), fill = Corr_value)) +
  geom_tile() +
  scale_fill_gradient2(low = "dodgerblue4", mid = "white", high = "green4", midpoint = 0) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(fill = "Correlation") +
  coord_fixed()
ggsave("Correlations/Heatmaps/fungal_root_genus_heatmap_v2.png", width = 5.5, height = 3.5, dpi = 300)


#### Fungal soil #####

# Get correlation value (R) and P value (P)
fun.soil.r <- rcorr(as.matrix(fun.soil.sig[2:ncol(fun.soil.sig)]))$r

# subset for just the relationships we want
fun.soil.r.format <- fun.soil.r %>% 
  as.data.frame(fun.soil.r)  %>%
  dplyr::select(1:(length(fun.soil.gen.names)-1)) %>%
  slice(-(c(1:(length(fun.soil.gen.names)-1), length(fun.soil.r[,1]))))
rownames(fun.soil.r.format) <- sigvar.names[1:13]

# pivot df longer
fun.soil_long <- fun.soil.r.format %>%
  rownames_to_column("Variable") %>%
  pivot_longer(cols = -Variable, names_to = "Genus", values_to = "Corr_value")

# make heatmap
fun.soil.Yorder <- c("Phialocephala", "Thelephora", "Athelopsis", "Diplodia", "Ascocorticium", "Tomentella", "Leohumicola", "Clavulinopsis", "Pseudocamarosporium", "Flavocillium")
ggplot(fun.soil_long, aes(x = factor(Variable, X.order), y = factor(Genus, fun.soil.Yorder), fill = Corr_value)) +
  geom_tile() +
  scale_fill_gradient2(low = "dodgerblue4", mid = "white", high = "green4", midpoint = 0) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(fill = "Correlation") +
  coord_fixed()
ggsave("Correlations/Heatmaps/fungal_soil_genus_heatmap_v2.png", width = 5.5, height = 3.5, dpi = 300)




