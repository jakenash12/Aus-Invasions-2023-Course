# Differentially abundant genus heatmap
# 1/6/24 CV

library(DescTools)
library(tidyverse)
library(tidyr)
library(ggplot2)

# read in full data
fulldat <- read.csv("Aus-Invasions-2023-Course/Merged_data/Aus23_allData_19Nov24.csv", row.names=1)

# read in genus-level pooled OTU tables
otu.genus.16S <- read.csv("Aus-Invasions-2023-Course/Aus23_16S_Metabarcoding/OTUTables/OtuMat16S_rare_byGenus_pooled.csv", row.names=1)
otu.genus.its <- read.csv("Aus-Invasions-2023-Course/Aus23_ITS_Metabarcoding/OTUTables/OtuMatITS_rare_byGenus_pooled.csv", row.names=1)
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

# pull out significant environmental variables and differentially abundant genera
sigvars <- c("treeName", "Tween_40_BiologDay5", "L.Serine_BiologDay5", "Glycyl.L.Glutamic_Acid_BiologDay5",
             "i.Erythritol_BiologDay5", "Putrescine_BiologDay5", "L.Asparagine_BiologDay5", 
             "pine_litter_prop", "Euc_leafLitter_percP", "Euc_leafLitter_totP",              
             "Litter_OLayer_ergosterol", "Litter_avg_ergosterol","euc_litter_prop", "Path_abund_root", "ergosterol","soil_moisture",
             "Bac_Shannon_soil")
sigdat <- fulldat[sigvars]

# merge each OTU table with environmental data
bac.root.sig <- merge(bac.root.gen, sigdat, by = "treeName")
fun.soil.sig <- merge(fun.soil.gen, sigdat, by = "treeName")

##################### Run correlations ##################### 
library(Hmisc)

#### Bacterial roots #####

# Get correlation value (R) and P value (P)
bac.root.r <- rcorr(as.matrix(bac.root.sig[2:ncol(bac.root.sig)]))$r
bac.root.p <- rcorr(as.matrix(bac.root.sig[2:ncol(bac.root.sig)]))$P

# subset for just the relationships we want
bac.root.r.format <- bac.root.r %>% 
  as.data.frame(bac.root.r)  %>%
  select(1:(length(bac.root.gen.names)-1)) %>% 
  slice(-(1:(length(bac.root.gen.names)-1)))
rownames(bac.root.r.format) <- sigvars[2:17]

bac.root.p.format <- bac.root.p %>% 
  as.data.frame(bac.root.p)  %>%
  select(1:(length(bac.root.gen.names)-1)) %>% 
  slice(-(1:(length(bac.root.gen.names)-1)))
rownames(bac.root.p.format) <- sigvars[2:17]

# pivot df longer
bac.root_long <- bac.root.r.format %>%
  rownames_to_column("Variable") %>%
  pivot_longer(cols = -Variable, names_to = "Genus", values_to = "Corr_value")

# add p-values
bac.rootP_long <- bac.root.p.format %>%
  rownames_to_column("Variable") %>%
  pivot_longer(cols = -Variable, names_to = "Genus", values_to = "P_value")
# generate stars to indicate significance
bac.rootP_long$P_star <- " "
for(i in 1:nrow(bac.rootP_long)){
  val <- bac.rootP_long$P_value[i]
  if(val < 0.001){
    bac.rootP_long$P_star[i] <- "***"
  } else if(val < 0.01){
    bac.rootP_long$P_star[i] <- "**"
  } else if(val < 0.05){
    bac.rootP_long$P_star[i] <- "*"
  }
}

bac.root_long$P_star <- bac.rootP_long$P_star

# make heatmap
ggplot(bac.root_long, aes(x = Genus, y = Variable, fill = Corr_value)) +
  geom_tile() +
  scale_fill_gradient2(low = "dodgerblue4", mid = "white", high = "green4", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(fill = "Correlation") +
  coord_fixed() + 
  geom_text(aes(label = P_star), size = 7)

#### Fungal soil #####

# Get correlation value (R) and P value (P)
fun.soil.r <- rcorr(as.matrix(fun.soil.sig[2:ncol(fun.soil.sig)]))$r
fun.soil.p <- rcorr(as.matrix(fun.soil.sig[2:ncol(fun.soil.sig)]))$P

# subset for just the relationships we want
fun.soil.r.format <- fun.soil.r %>% 
  as.data.frame(fun.soil.r)  %>%
  select(1:(length(fun.soil.gen.names)-1)) %>%
  slice(-(1:length(fun.soil.gen.names)-1))
rownames(fun.soil.r.format) <- sigvars[2:17]

fun.soil.p.format <- fun.soil.p %>% 
  as.data.frame(fun.soil.p)  %>%
  select(1:(length(fun.soil.gen.names)-1)) %>%
  slice(-(1:length(fun.soil.gen.names)-1))
rownames(fun.soil.p.format) <- sigvars[2:17]

# pivot df longer
fun.soil_long <- fun.soil.r.format %>%
  rownames_to_column("Variable") %>%
  pivot_longer(cols = -Variable, names_to = "Genus", values_to = "Corr_value")

# add p-values
fun.soilP_long <- fun.soil.p.format %>%
  rownames_to_column("Variable") %>%
  pivot_longer(cols = -Variable, names_to = "Genus", values_to = "P_value")
# generate stars to indicate significance
fun.soilP_long$P_star <- " "
for(i in 1:nrow(fun.soilP_long)){
  val <- fun.soilP_long$P_value[i]
  if(val < 0.001){
    fun.soilP_long$P_star[i] <- "***"
  } else if(val < 0.01){
    fun.soilP_long$P_star[i] <- "**"
  } else if(val < 0.05){
    fun.soilP_long$P_star[i] <- "*"
  }
}

fun.soil_long$P_star <- fun.soilP_long$P_star

# make heatmap
ggplot(fun.soil_long, aes(x = Genus, y = Variable, fill = Corr_value)) +
  geom_tile() +
  scale_fill_gradient2(low = "dodgerblue4", mid = "white", high = "green4", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(fill = "Correlation") +
  coord_fixed() + 
  geom_text(aes(label = P_star), size = 7)

