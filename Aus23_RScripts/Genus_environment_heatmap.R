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

colnames(soil.16s)[4] <- "11-24"

# subset out for the genera we want
bac.root.gen.names <- c("treeName", "Lechevalieria", "Nocardia", "Acidicaldus", "Nevskia", "Occallatibacter", "SBR1031", "Novosphingobium", "Collimonas", "A4b", "Bryocella")
bac.soil.gen.names <- c("treeName", "Gryllotalpicola", "Limnobacter", "Acidisoma", "Dinghuibacter", "Acidiphilium", "11-24", "SM1A02", "Subgroup_5", "A4b", "Mesorhizobium")
# "uncultured_24" "11-24", uncultured_40, uncultured_23
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
sigvar.order <-c("Tween_40_BiologDay5", "L.Serine_BiologDay5", 
    "Glycyl.L.Glutamic_Acid_BiologDay5",
    "i.Erythritol_BiologDay5", "Putrescine_BiologDay5", "L.Asparagine_BiologDay5", "soil_moisture", "ergosterol", 
    "Litter_OLayer_ergosterol", "Litter_avg_ergosterol",
    "pine_litter_prop", "euc_litter_prop", "Euc_leafLitter_percP")
sigdat <- fulldat[sigvars]

# merge each OTU table with environmental data
bac.root.sig <- merge(bac.root.gen, sigdat, by = "treeName")
bac.soil.sig <- merge(bac.soil.gen, sigdat, by = "treeName")
fun.root.sig <- merge(fun.root.gen, sigdat, by = "treeName")
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
  slice(-(c(1:(length(bac.root.gen.names)-1), length(bac.root.r[,1]))))
rownames(bac.root.r.format) <- sigvars[2:14]

bac.root.p.format <- bac.root.p %>% 
  as.data.frame(bac.root.p)  %>%
  select(1:(length(bac.root.gen.names)-1)) %>% 
  slice(-(c(1:(length(bac.root.gen.names)-1), length(bac.root.r[,1]))))
rownames(bac.root.p.format) <- sigvars[2:14]

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
bac.root.Yorder <- c("Lechevalieria", "Nocardia", "Acidicaldus", "Nevskia", "Occallatibacter", "SBR1031", "Novosphingobium", "Collimonas", "A4b", "Bryocella")
ggplot(bac.root_long, aes(x = factor(Variable, sigvar.order), y = factor(Genus, bac.root.Yorder), fill = Corr_value)) +
  geom_tile() +
  scale_fill_gradient2(low = "dodgerblue4", mid = "white", high = "green4", midpoint = 0) +
  labs(x = "Environmental variable", y = "Bacterial genus") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(fill = "Correlation") +
  coord_fixed() + 
  geom_text(aes(label = P_star), size = 5)
ggsave("/Users/crvietorisz/Desktop/GitHub/Aus-Invasions-2023-Course/Correlations/Heatmaps/bacteria_root_genus_heatmap.png", width = 6, heigh = 5, dpi = 300)

#### Bacterial soil #####

# Get correlation value (R) and P value (P)
bac.soil.r <- rcorr(as.matrix(bac.soil.sig[2:ncol(bac.soil.sig)]))$r
bac.soil.p <- rcorr(as.matrix(bac.soil.sig[2:ncol(bac.soil.sig)]))$P

# subset for just the relationships we want
bac.soil.r.format <- bac.soil.r %>% 
  as.data.frame(bac.soil.r)  %>%
  select(c(1:(length(bac.soil.gen.names)-1), length(bac.soil.r[,1]))) %>% 
  slice(-(c(1:(length(bac.soil.gen.names)-1), length(bac.soil.r[,1]))))
rownames(bac.soil.r.format) <- sigvars[2:14]

bac.soil.p.format <- bac.soil.p %>% 
  as.data.frame(bac.soil.p)  %>%
  select(c(1:(length(bac.soil.gen.names)-1), length(bac.soil.r[,1]))) %>% 
  slice(-(c(1:(length(bac.soil.gen.names)-1), length(bac.soil.r[,1]))))
rownames(bac.soil.p.format) <- sigvars[2:14]

# pivot df longer
bac.soil_long <- bac.soil.r.format %>%
  rownames_to_column("Variable") %>%
  pivot_longer(cols = -Variable, names_to = "Genus", values_to = "Corr_value")

# add p-values
bac.soilP_long <- bac.soil.p.format %>%
  rownames_to_column("Variable") %>%
  pivot_longer(cols = -Variable, names_to = "Genus", values_to = "P_value")
# generate stars to indicate significance
bac.soilP_long$P_star <- " "
for(i in 1:nrow(bac.soilP_long)){
  val <- bac.soilP_long$P_value[i]
  if(val < 0.001){
    bac.soilP_long$P_star[i] <- "***"
  } else if(val < 0.01){
    bac.soilP_long$P_star[i] <- "**"
  } else if(val < 0.05){
    bac.soilP_long$P_star[i] <- "*"
  }
}

bac.soil_long$P_star <- bac.soilP_long$P_star

# make heatmap
bac.soil.Yorder <- c("Gryllotalpicola", "Limnobacter", "Acidisoma", "Dinghuibacter", "SM1A02", "Subgroup_5", "Bac_Shannon_soil")
ggplot(bac.soil_long, aes(x = factor(Variable, sigvar.order), y = factor(Genus, bac.soil.Yorder), fill = Corr_value)) +
  geom_tile() +
  scale_fill_gradient2(low = "dodgerblue4", mid = "white", high = "green4", midpoint = 0) +
  labs(x = "Environmental variable", y = "Bacterial genus or diversity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(fill = "Correlation") +
  coord_fixed() + 
  geom_text(aes(label = P_star), size = 5)
ggsave("/Users/crvietorisz/Desktop/GitHub/Aus-Invasions-2023-Course/Correlations/Heatmaps/bacteria_soil_genus_heatmap.png", width = 6, heigh = 5, dpi = 300)


#### Fungal roots #####

# Get correlation value (R) and P value (P)
fun.root.r <- rcorr(as.matrix(fun.root.sig[2:ncol(fun.root.sig)]))$r
fun.root.p <- rcorr(as.matrix(fun.root.sig[2:ncol(fun.root.sig)]))$P

# subset for just the relationships we want
fun.root.r.format <- fun.root.r %>% 
  as.data.frame(fun.root.r)  %>%
  select(1:(length(fun.root.gen.names)-1)) %>%
  slice(-(c(1:(length(fun.root.gen.names)-1), length(fun.root.r[,1]))))
rownames(fun.root.r.format) <- sigvars[2:14]

fun.root.p.format <- fun.root.p %>% 
  as.data.frame(fun.root.p)  %>%
  select(1:(length(fun.root.gen.names)-1)) %>%
  slice(-(c(1:(length(bac.root.gen.names)-1), length(bac.root.r[,1]))))
rownames(fun.root.p.format) <- sigvars[2:14]

# pivot df longer
fun.root_long <- fun.root.r.format %>%
  rownames_to_column("Variable") %>%
  pivot_longer(cols = -Variable, names_to = "Genus", values_to = "Corr_value")

# add p-values
fun.rootP_long <- fun.root.p.format %>%
  rownames_to_column("Variable") %>%
  pivot_longer(cols = -Variable, names_to = "Genus", values_to = "P_value")
# generate stars to indicate significance
fun.rootP_long$P_star <- " "
for(i in 1:nrow(fun.rootP_long)){
  val <- fun.rootP_long$P_value[i]
  if(val < 0.001){
    fun.rootP_long$P_star[i] <- "***"
  } else if(val < 0.01){
    fun.rootP_long$P_star[i] <- "**"
  } else if(val < 0.05){
    fun.rootP_long$P_star[i] <- "*"
  }
}

fun.root_long$P_star <- fun.rootP_long$P_star

# make heatmap
fun.root.Yorder <- c("Phialocephala", "Talaromyces", "Devriesia", "Aspergillus", "Penicillium", "Auricularia", "Pezoloma", "Cryptosporiopsis", "Pezicula", "Mollisia")
ggplot(fun.root_long, aes(x = factor(Variable, sigvar.order), y = factor(Genus, fun.root.Yorder), fill = Corr_value)) +
  geom_tile() +
  scale_fill_gradient2(low = "dodgerblue4", mid = "white", high = "green4", midpoint = 0) +
  labs(x = "Environmental variable", y = "Fungal genus") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(fill = "Correlation") +
  coord_fixed() + 
  geom_text(aes(label = P_star), size = 5)
ggsave("/Users/crvietorisz/Desktop/GitHub/Aus-Invasions-2023-Course/Correlations/Heatmaps/fungal_root_genus_heatmap.png", width = 6, heigh = 5, dpi = 300)


#### Fungal soil #####

# Get correlation value (R) and P value (P)
fun.soil.r <- rcorr(as.matrix(fun.soil.sig[2:ncol(fun.soil.sig)]))$r
fun.soil.p <- rcorr(as.matrix(fun.soil.sig[2:ncol(fun.soil.sig)]))$P

# subset for just the relationships we want
fun.soil.r.format <- fun.soil.r %>% 
  as.data.frame(fun.soil.r)  %>%
  select(1:(length(fun.soil.gen.names)-1)) %>%
  slice(-(c(1:(length(fun.soil.gen.names)-1), length(fun.soil.r[,1]))))
rownames(fun.soil.r.format) <- sigvars[2:14]

fun.soil.p.format <- fun.soil.p %>% 
  as.data.frame(fun.soil.p)  %>%
  select(1:(length(fun.soil.gen.names)-1)) %>%
  slice(-(c(1:(length(bac.root.gen.names)-1), length(bac.root.r[,1]))))
rownames(fun.soil.p.format) <- sigvars[2:14]

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
fun.soil.Yorder <- c("Phialocephala", "Thelephora", "Athelopsis", "Diplodia", "Ascocorticium", "Tomentella", "Leohumicola", "Clavulinopsis", "Pseudocamarosporium", "Flavocillium")
ggplot(fun.soil_long, aes(x = factor(Variable, sigvar.order), y = factor(Genus, fun.soil.Yorder), fill = Corr_value)) +
  geom_tile() +
  scale_fill_gradient2(low = "dodgerblue4", mid = "white", high = "green4", midpoint = 0) +
  labs(x = "Environmental variable", y = "Fungal genus") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(fill = "Correlation") +
  coord_fixed() + 
  geom_text(aes(label = P_star), size = 7)
ggsave("/Users/crvietorisz/Desktop/GitHub/Aus-Invasions-2023-Course/Correlations/Heatmaps/fungal_soil_genus_heatmap.png", width = 6, heigh = 5, dpi = 300)




