# Create a stacked barplot of 18S taxa relative abundances by tree species
# C. Vietorisz 12/29/24

library(DescTools)
library(tidyverse)
library(phyloseq)


# read in otu tables
amf.otu <- read.csv("Aus-Invasions-2023-Course/Aus23_16S_Metabarcoding/OTUTables/OtuMat16S_rel.csv", row.names=1)
amf.soil.otu <- subset(amf.otu, rownames(amf.otu) %like% "%Soil%")
amf.root.otu <- subset(amf.otu, rownames(amf.otu) %like% "%Root%")

# read in taxonomy list
amf.tax <- read_tsv("Aus-Invasions-2023-Course/Aus23_16S_Metabarcoding/Aus23_16S_SilvaTaxonomy_16S.tsv")[1:2]
amf.tax$Taxon <- str_replace_all(amf.tax$Taxon, "\\w__", "") # remove 
amf.tax.split <- amf.tax %>% 
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
  column_to_rownames(var = "Feature ID")


######## SOIL #########
# put amfk into phyloseq object
ps.amf.soil <- phyloseq(otu_table(as.matrix(amf.soil.otu), taxa_are_rows=FALSE), 
                        tax_table(as.matrix(amf.tax.split)))
# melt OTU table so each taxon is a row
psz.amf.soil <- psmelt(ps.amf.soil)
# add tree species
psz.amf.soil$Tree <- ifelse(grepl("E", psz.amf.soil$Sample), "Eucalypt",
                            ifelse(grepl("P", psz.amf.soil$Sample), "Pine", NA))
write.csv(psz.amf.soil, "Aus-Invasions-2023-Course/Aus23_16S_Metabarcoding/OTUTables/OTU_16S_soil_rel_melted.csv")


#Select top 20 most abundant families
amf.soil.20fam.names <- psz.amf.soil %>%
  filter(!is.na(Family) & Family != "NA" & Family != " uncultured") %>%
  group_by(Family) %>%
  summarize(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 20) %>%
  pull(Family) 
amf.soil.20fam.names
amf.soil.20fam <- subset(psz.amf.soil, psz.amf.soil$Family %like% amf.soil.20fam.names)

######## ROOTS #########
# put amfk into phyloseq object
ps.amf.root <- phyloseq(otu_table(as.matrix(amf.root.otu), taxa_are_rows=FALSE), 
                        tax_table(as.matrix(amf.tax.split)))
# melt OTU table so each taxon is a row
psz.amf.root <- psmelt(ps.amf.root)
# add tree species
psz.amf.root$Tree <- ifelse(grepl("E", psz.amf.root$Sample), "Eucalypt",
                            ifelse(grepl("P", psz.amf.root$Sample), "Pine", NA))
write.csv(psz.amf.root, "Aus-Invasions-2023-Course/Aus23_16S_Metabarcoding/OTUTables/OTU_16S_root_rel_melted.csv")

#Select top 20 most abundant genera
amf.root.20fam.names <- psz.amf.root %>%
  filter(!is.na(Family) & Family != "NA" & Family != " uncultured" & Family != " Unknown_Family") %>%
  group_by(Family) %>%
  summarize(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 20) %>%
  pull(Family) 
amf.root.20fam.names
amf.root.20fam <- subset(psz.amf.root, psz.amf.root$Family %like% amf.root.20fam.names)


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

amf.soil.plot <- ggplot(amf.soil.20fam, aes(y=Abundance, x=Tree, fill = Family, color = Family)) + 
  geom_bar(position="fill", stat="identity")+
  labs(x = "Tree species", y = "Relative abundance")+
  scale_fill_manual(values = colors)+ 
  scale_color_manual(values = colors)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.amfkground = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13), axis.text=element_text(size=12))
amf.soil.plot
ggsave("16S_soil_stackedBarplot_top20genera.png", path = "Aus-Invasions-2023-Course/Aus23_16S_Metabarcoding/StackedBarplots", width=8, height=10, dpi=300)

amf.root.plot <- ggplot(amf.root.20fam, aes(y=Abundance, x=Tree, fill = Family, color = Family)) + 
  geom_bar(position="fill", stat="identity")+
  labs(x = "Tree species", y = "Relative abundance")+
  scale_fill_manual(values = colors)+
  scale_color_manual(values = colors)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.amfkground = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13), axis.text=element_text(size=12))
amf.root.plot
ggsave("16S_root_stackedBarplot_top20genera.png", path = "Aus-Invasions-2023-Course/Aus23_16S_Metabarcoding/StackedBarplots", width=8, height=10, dpi=300)







