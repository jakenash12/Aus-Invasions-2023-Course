# Create a stacked barplot of 18S taxa relative abundances by tree species
# C. Vietorisz 12/29/24

library(DescTools)
library(tidyverse)
library(phyloseq)


# read in otu tables
amf.otu <- read.csv("Aus-Invasions-2023-Course/Aus23_18S_Metabarcoding/OTUTables/OtuMat18S_rel.csv", row.names=1)
amf.soil.otu <- subset(amf.otu, rownames(amf.otu) %like% "%Soil%")
amf.root.otu <- subset(amf.otu, rownames(amf.otu) %like% "%Root%")

# read in taxonomy list
amf.tax <- read.csv("Aus-Invasions-2023-Course/Aus23_18S_Metabarcoding/AMF_Taxonomy_SplitFormat.csv")
rownames(amf.tax) <- amf.tax[,1]
amf.tax <- amf.tax[2:8]

######## SOIL #########
# put back into phyloseq object
ps.amf.soil <- phyloseq(otu_table(as.matrix(amf.soil.otu), taxa_are_rows=FALSE), 
                        tax_table(as.matrix(amf.tax)))
# melt OTU table so each taxon is a row
psz.amf.soil <- psmelt(ps.amf.soil)
# add tree species
psz.amf.soil$Tree <- ifelse(grepl("E", psz.amf.soil$Sample), "Eucalypt",
                            ifelse(grepl("P", psz.amf.soil$Sample), "Pine", NA))
write.csv(psz.amf.soil, "Aus-Invasions-2023-Course/Aus23_18S_Metabarcoding/OTUTables/OTU_18S_soil_rel_melted.csv")


######## ROOTS #########
# put back into phyloseq object
ps.amf.root <- phyloseq(otu_table(as.matrix(amf.root.otu), taxa_are_rows=FALSE), 
                        tax_table(as.matrix(amf.tax)))
# melt OTU table so each taxon is a row
psz.amf.root <- psmelt(ps.amf.root)
# add tree species
psz.amf.root$Tree <- ifelse(grepl("E", psz.amf.root$Sample), "Eucalypt",
                            ifelse(grepl("P", psz.amf.root$Sample), "Pine", NA))
write.csv(psz.amf.root, "Aus-Invasions-2023-Course/Aus23_18S_Metabarcoding/OTUTables/OTU_18S_root_rel_melted.csv")

############ SOIL AND ROOTS ############ 

# put back into phyloseq object
ps.amf <- phyloseq(otu_table(as.matrix(amf.otu), taxa_are_rows=FALSE), 
                   tax_table(as.matrix(amf.tax)))
# melt OTU table so each taxon is a row
psz.amf <- psmelt(ps.amf)
# add tree species
psz.amf$SampleCategory <- ifelse(grepl("E", psz.amf$Sample),
                                 ifelse(grepl("Soil", psz.amf$Sample), "Eucalypt Soil",
                                        ifelse(grepl("Root", psz.amf$Sample), "Eucalypt Root", NA)),
                                 ifelse(grepl("P", psz.amf$Sample),
                                        ifelse(grepl("Soil", psz.amf$Sample), "Pine Soil",
                                               ifelse(grepl("Root", psz.amf$Sample), "Pine Root", NA)),
                                        NA))
write.csv(psz.amf, "Aus-Invasions-2023-Course/Aus23_18S_Metabarcoding/OTUTables/OTU_18S_SoilRoot_rel_melted.csv")


#################################### MAKE PLOTS #################################### 

# get colors
colors5 <- c(
  "#1f77b4", # Steel Blue
  "#2ca02c", # Forest Green
  "#ff7f0e", # Dark Orange
  "#9edae5", # Pale Blue
  "#9467bd" # Medium Purple
)
######

amf.soil.plot <- ggplot(psz.amf.soil, aes(y=Abundance, x=Tree, fill = Genus, color = Genus)) + 
  geom_bar(position="fill", stat="identity")+
  labs(x = "Tree species", y = "Relative abundance")+
  scale_fill_manual(values = colors5)+ 
  scale_color_manual(values = colors5)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13), axis.text=element_text(size=12))
amf.soil.plot
ggsave("18S_soil_stackedBarplot_genus.png", path = "Aus-Invasions-2023-Course/Aus23_18S_Metabarcoding/StackedBarplots", width=8, height=10, dpi=300)

amf.root.plot <- ggplot(psz.amf.root, aes(y=Abundance, x=Tree, fill = Genus, color = Genus)) + 
  geom_bar(position="fill", stat="identity")+
  labs(x = "Tree species", y = "Relative abundance")+
  scale_fill_manual(values = colors5)+
  scale_color_manual(values = colors5)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13), axis.text=element_text(size=12))
amf.root.plot
ggsave("18S_root_stackedBarplot_genus.png", path = "Aus-Invasions-2023-Course/Aus23_18S_Metabarcoding/StackedBarplots", width=8, height=10, dpi=300)


cat.ord <- c("Eucalypt Root", "Pine Root", "Eucalypt Soil", "Pine Soil")
amf.all.plot <- ggplot(psz.amf, aes(y=Abundance, x=factor(SampleCategory, cat.ord), fill = Genus, color = Genus)) + 
  geom_bar(position="fill", stat="identity")+
  labs(x = "Sample type", y = "Relative abundance")+
  scale_fill_manual(values = colors5)+ 
  scale_color_manual(values = colors5)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13), axis.text=element_text(size=12))
amf.all.plot
ggsave("18S_SoilRoot_stackedBarplot_genus.png", path = "Aus-Invasions-2023-Course/Aus23_18S_Metabarcoding/StackedBarplots", width=7, height=6, dpi=300)





