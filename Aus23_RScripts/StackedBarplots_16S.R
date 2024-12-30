# Create a stacked barplot of 16S taxa relative abundances by tree species
# C. Vietorisz 12/29/24

library(DescTools)
library(tidyverse)
library(phyloseq)


# read in otu tables
bac.otu <- read.csv("Aus-Invasions-2023-Course/Aus23_16S_Metabarcoding/OTUTables/OtuMat16S_rel.csv", row.names=1)
bac.soil.otu <- subset(bac.otu, rownames(bac.otu) %like% "%Soil%")
bac.root.otu <- subset(bac.otu, rownames(bac.otu) %like% "%Root%")

# read in taxonomy list
bac.tax <- read_tsv("Aus-Invasions-2023-Course/Aus23_16S_Metabarcoding/Aus23_16S_SilvaTaxonomy_16S.tsv")[1:2]
bac.tax$Taxon <- str_replace_all(bac.tax$Taxon, "\\w__", "") # remove 
bac.tax.split <- bac.tax %>% 
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
  column_to_rownames(var = "Feature ID")
  

######## SOIL #########
# put back into phyloseq object
ps.bac.soil <- phyloseq(otu_table(as.matrix(bac.soil.otu), taxa_are_rows=FALSE), 
                        tax_table(as.matrix(bac.tax.split)))
# melt OTU table so each taxon is a row
psz.bac.soil <- psmelt(ps.bac.soil)
# add tree species
psz.bac.soil$Tree <- ifelse(grepl("E", psz.bac.soil$Sample), "Eucalypt",
                            ifelse(grepl("P", psz.bac.soil$Sample), "Pine", NA))
write.csv(psz.bac.soil, "Aus-Invasions-2023-Course/Aus23_16S_Metabarcoding/OTUTables/OTU_16S_soil_rel_melted.csv")


#Select top 20 most abundant families
bac.soil.20fam.names <- psz.bac.soil %>%
  filter(!is.na(Family) & Family != "NA" & Family != " uncultured") %>%
  group_by(Family) %>%
  summarize(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 20) %>%
  pull(Family) 
bac.soil.20fam.names
bac.soil.20fam <- subset(psz.bac.soil, psz.bac.soil$Family %like% bac.soil.20fam.names)

######## ROOTS #########
# put back into phyloseq object
ps.bac.root <- phyloseq(otu_table(as.matrix(bac.root.otu), taxa_are_rows=FALSE), 
                        tax_table(as.matrix(bac.tax.split)))
# melt OTU table so each taxon is a row
psz.bac.root <- psmelt(ps.bac.root)
# add tree species
psz.bac.root$Tree <- ifelse(grepl("E", psz.bac.root$Sample), "Eucalypt",
                            ifelse(grepl("P", psz.bac.root$Sample), "Pine", NA))
write.csv(psz.bac.root, "Aus-Invasions-2023-Course/Aus23_16S_Metabarcoding/OTUTables/OTU_16S_root_rel_melted.csv")

#Select top 20 most abundant genera
bac.root.20fam.names <- psz.bac.root %>%
  filter(!is.na(Family) & Family != "NA" & Family != " uncultured" & Family != " Unknown_Family") %>%
  group_by(Family) %>%
  summarize(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 20) %>%
  pull(Family) 
bac.root.20fam.names
bac.root.20fam <- subset(psz.bac.root, psz.bac.root$Family %like% bac.root.20fam.names)

############ SOIL AND ROOTS ############ 

# put back into phyloseq object
ps.bac <- phyloseq(otu_table(as.matrix(bac.otu), taxa_are_rows=FALSE), 
                   tax_table(as.matrix(bac.tax.split)))

# melt OTU table so each taxon is a row
psz.bac <- psmelt(ps.bac)
# add tree species
psz.bac$SampleCategory <- ifelse(grepl("E", psz.bac$Sample),
                                 ifelse(grepl("Soil", psz.bac$Sample), "Eucalypt Soil",
                                        ifelse(grepl("Root", psz.bac$Sample), "Eucalypt Root", NA)),
                                 ifelse(grepl("P", psz.bac$Sample),
                                        ifelse(grepl("Soil", psz.bac$Sample), "Pine Soil",
                                               ifelse(grepl("Root", psz.bac$Sample), "Pine Root", NA)),
                                        NA))
write.csv(psz.bac, "Aus-Invasions-2023-Course/Aus23_16S_Metabarcoding/OTUTables/OTU_16S_SoilRoot_rel_melted.csv")


bac.30fam.names <- psz.bac %>%
  filter(!is.na(Family) & Family != "NA" & Family != " uncultured" & Family != " Unknown_Family") %>%
  group_by(Family) %>%
  summarize(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 30) %>%
  pull(Family) 
bac.30fam.names
bac.30fam <- subset(psz.bac, psz.bac$Family %like% bac.30fam.names & !is.na(SampleCategory))




#################################### MAKE PLOTS #################################### 

colors <- c(
  "#1f77b4", # Steel Blue
  "#ff7f0e", # Dark Orange
  "#2ca02c", # Forest Green
  "#d62728", # Brick Red
  "#9467bd", # Medium Purple
  "#8c564b", # Brown
  "#e377c2", # Pink
  "#7f7f7f", # Gray
  "#bcbd22", # Olive
  "#17becf", # Light Blue
  "#ff9896", # Light Salmon
  "#98df8a", # Light Green
  "#ffbb78", # Light Orange
  "#c5b0d5", # Light Purple
  "#c49c94", # Light Brown
  "#f7b6d2", # Light Pink
  "#dbdb8d", # Khaki
  "#9edae5", # Pale Blue
  "#393b79", # Dark Blue
  "#7b4173"  # Dark Purple
)

bac.soil.plot <- ggplot(bac.soil.20fam, aes(y=Abundance, x=Tree, fill = Family, color = Family)) + 
  geom_bar(position="fill", stat="identity")+
  labs(x = "Tree species", y = "Relative abundance")+
  scale_fill_manual(values = colors)+ 
  scale_color_manual(values = colors)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13), axis.text=element_text(size=12))
bac.soil.plot
ggsave("16S_soil_stackedBarplot_top20family.png", path = "Aus-Invasions-2023-Course/Aus23_16S_Metabarcoding/StackedBarplots", width=8, height=10, dpi=300)

bac.root.plot <- ggplot(bac.root.20fam, aes(y=Abundance, x=Tree, fill = Family, color = Family)) + 
  geom_bar(position="fill", stat="identity")+
  labs(x = "Tree species", y = "Relative abundance")+
  scale_fill_manual(values = colors)+
  scale_color_manual(values = colors)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13), axis.text=element_text(size=12))
bac.root.plot
ggsave("16S_root_stackedBarplot_top20family.png", path = "Aus-Invasions-2023-Course/Aus23_16S_Metabarcoding/StackedBarplots", width=8, height=10, dpi=300)

cat.ord <- c("Eucalypt Root", "Pine Root", "Eucalypt Soil", "Pine Soil")
bac.all.plot <- ggplot(bac.30fam, aes(y=Abundance, x=factor(SampleCategory, cat.ord), fill = Family, color = Family)) + 
  geom_bar(position="fill", stat="identity")+
  labs(x = "Sample type", y = "Relative abundance")+
  scale_fill_manual(values = colors30)+ 
  scale_color_manual(values = colors30)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13), axis.text=element_text(size=12))
bac.all.plot
ggsave("16S_SoilRoot_stackedBarplot_top30family.png", path = "Aus-Invasions-2023-Course/Aus23_16S_Metabarcoding/StackedBarplots", width=9, height=8, dpi=300)






