# Create a stacked barplot of taxa relative abundances by tree species
# C. Vietorisz 12/29/24

library(DescTools)
library(tidyverse)
library(phyloseq)

setwd("/Users/moniquegagnon/Desktop/GitHub/")

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

# put back into phyloseq object
ps.its.soil <- phyloseq(otu_table(as.matrix(its.soil.otu), taxa_are_rows=FALSE), 
               tax_table(as.matrix(its.tax.split)))

# melt OTU table so each taxon is a row
psz.its.soil <- psmelt(ps.its.soil)
# add tree species
psz.its.soil$Tree <- ifelse(grepl("E", psz.its.soil$Sample), "Eucalypt",
                               ifelse(grepl("P", psz.its.soil$Sample), "Pine", NA))
write.csv(psz.its.soil, "Aus-Invasions-2023-Course/Aus23_ITS_Metabarcoding/OTUTables/OTU_ITS_rel_melted.csv")

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
  filter(!is.na(Genus) & Genus != "NA") %>%
  group_by(Genus) %>%
  summarize(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 20) %>%
  pull(Genus) 
its.soil.20gen <- subset(psz.its.soil, psz.its.soil$Genus %like% its.soil.20gen.names)

#################################### MAKE PLOTS #################################### 

its.soil.plot <- ggplot(its.soil.20gen, aes(y=Abundance, x=Tree, fill = Genus)) + 
  geom_bar(position="fill", stat="identity")+
  labs(x = "Tree species", y = "Relative abundance")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13), axis.text=element_text(size=12))
its.soil.plot


ggsave("EMF_SAP_guild_abund_by_forest_v2.png", path = "/Users/moniquegagnon/Desktop/BU/PhD/White pine/Figures/Soil microbes/Taxa abunds", width=9, height=5, dpi=300)



