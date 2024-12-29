# Create a stacked barplot of ITS taxa relative abundances by tree species
# C. Vietorisz 12/29/24

library(DescTools)
library(tidyverse)
library(phyloseq)

# read in otu tables
its.otu <- read.csv("Aus-Invasions-2023-Course/Aus23_ITS_Metabarcoding/OTUTables/OtuMatITS_rel.csv", row.names=1)
its.soil.otu <- subset(its.otu, rownames(its.otu) %like% "%Soil%")
its.root.otu <- subset(its.otu, rownames(its.otu) %like% "%Root%")

# read in taxonomy list
its.tax <- read_tsv("Aus-Invasions-2023-Course/Aus23_ITS_Metabarcoding/ITS2_Dada2_repseqs97_taxonomy_edited_FunT_FunG.tsv")[c(1,6,7)]
its.tax.split <- its.tax %>% 
  separate(Kingdom, into = c("Phylum", "Class", "Order", "Family", "Genus"), sep = ";") %>%
  column_to_rownames(var = "Feature ID") %>%
  mutate(Kingdom = "Fungi") %>%
  select(Kingdom, everything())

######## SOIL #########
# put back into phyloseq object
ps.its.soil <- phyloseq(otu_table(as.matrix(its.soil.otu), taxa_are_rows=FALSE), 
               tax_table(as.matrix(its.tax.split)))

# melt OTU table so each taxon is a row
psz.its.soil <- psmelt(ps.its.soil)
# add tree species
psz.its.soil$Tree <- ifelse(grepl("E", psz.its.soil$Sample), "Eucalypt",
                               ifelse(grepl("P", psz.its.soil$Sample), "Pine", NA))
write.csv(psz.its.soil, "Aus-Invasions-2023-Course/Aus23_ITS_Metabarcoding/OTUTables/OTU_ITS_soil_rel_melted.csv")

#Selecting the Top 100 most abundant taxa
its.soil.top100 <- names(sort(taxa_sums(ps.its.soil), decreasing=TRUE))[1:100]
ps100.its.soil <- prune_taxa(its.soil.top100, ps.its.soil)
# melt OTU table so each taxon is a row
psz100.its.soil <- psmelt(ps100.its.soil)
# add tree species
psz100.its.soil$Tree <- ifelse(grepl("E", psz100.its.soil$Sample), "Eucalypt",
                        ifelse(grepl("P", psz100.its.soil$Sample), "Pine", NA))

#Select top 20 most abundant genera
its.soil.20gen.names <- psz.its.soil %>%
  filter(!is.na(Genus) & Genus != "NA" & !grepl("Incertae_sedis", Genus)) %>%
  group_by(Genus) %>%
  summarize(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 20) %>%
  pull(Genus) 
its.soil.20gen <- subset(psz.its.soil, psz.its.soil$Genus %like% its.soil.20gen.names)

######## ROOTS #########
# put back into phyloseq object
ps.its.root <- phyloseq(otu_table(as.matrix(its.root.otu), taxa_are_rows=FALSE), 
                        tax_table(as.matrix(its.tax.split)))

# melt OTU table so each taxon is a row
psz.its.root <- psmelt(ps.its.root)
# add tree species
psz.its.root$Tree <- ifelse(grepl("E", psz.its.root$Sample), "Eucalypt",
                            ifelse(grepl("P", psz.its.root$Sample), "Pine", NA))
write.csv(psz.its.root, "Aus-Invasions-2023-Course/Aus23_ITS_Metabarcoding/OTUTables/OTU_ITS_root_rel_melted.csv")

#Select top 20 most abundant genera
its.root.20gen.names <- psz.its.root %>%
  filter(!is.na(Genus) & Genus != "NA" & !grepl("Incertae_sedis", Genus)) %>%
  group_by(Genus) %>%
  summarize(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 20) %>%
  pull(Genus) 
its.root.20gen <- subset(psz.its.root, psz.its.root$Genus %like% its.root.20gen.names)


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

its.soil.plot <- ggplot(its.soil.20gen, aes(y=Abundance, x=Tree, fill = Genus, color = Genus)) + 
  geom_bar(position="fill", stat="identity")+
  labs(x = "Tree species", y = "Relative abundance")+
  scale_fill_manual(values = colors)+
  scale_color_manual(values = colors)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13), axis.text=element_text(size=12))
its.soil.plot
ggsave("its_soil_stackedBarplot_top20genera.png", path = "Aus-Invasions-2023-Course/Aus23_ITS_Metabarcoding/StackedBarplots", width=8, height=10, dpi=300)

its.root.plot <- ggplot(its.root.20gen, aes(y=Abundance, x=Tree, fill = Genus, color = Genus)) + 
  geom_bar(position="fill", stat="identity")+
  labs(x = "Tree species", y = "Relative abundance")+
  scale_fill_manual(values = colors)+ 
  scale_color_manual(values = colors)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13), axis.text=element_text(size=12))
its.root.plot
ggsave("its_root_stackedBarplot_top20genera.png", path = "Aus-Invasions-2023-Course/Aus23_ITS_Metabarcoding/StackedBarplots", width=8, height=10, dpi=300)







