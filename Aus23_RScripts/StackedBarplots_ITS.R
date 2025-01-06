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
# append with info from Rytas on native vs introduced fungi
origin <- read.csv("Aus-Invasions-2023-Course/Aus23_ITS_Metabarcoding/ITS2_Dada2_repseqs97_RytasEdits.csv")[c(1,15)]
origin[origin ==""] <- NA
origin[origin =="introduced"] <- "Introduced"
colnames(origin)[1] <- "OTU"
its.tax.split$OTU <- rownames(its.tax.split)
its.tax.split <- merge(its.tax.split, origin, by = "OTU", all.x=TRUE)

############################################################################################################
# ALL FUNGI
############################################################################################################

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

############ SOIL AND ROOTS ############ 

# put back into phyloseq object
ps.its <- phyloseq(otu_table(as.matrix(its.otu), taxa_are_rows=FALSE), 
                        tax_table(as.matrix(its.tax.split)))

# melt OTU table so each taxon is a row
psz.its <- psmelt(ps.its)
# add tree species
psz.its$SampleCategory <- ifelse(grepl("E", psz.its$Sample),
                                      ifelse(grepl("Soil", psz.its$Sample), "Eucalypt Soil",
                                             ifelse(grepl("Root", psz.its$Sample), "Eucalypt Root", NA)),
                                      ifelse(grepl("P", psz.its$Sample),
                                             ifelse(grepl("Soil", psz.its$Sample), "Pine Soil",
                                                    ifelse(grepl("Root", psz.its$Sample), "Pine Root", NA)),
                                             NA))
write.csv(psz.its, "Aus-Invasions-2023-Course/Aus23_ITS_Metabarcoding/OTUTables/OTU_ITS_SoilRoot_rel_melted.csv")
psz.its <- read.csv("Aus-Invasions-2023-Course/Aus23_ITS_Metabarcoding/OTUTables/OTU_ITS_SoilRoot_rel_melted.csv", row.names=1)

# add native vs introduced column
psz.its <- merge(psz.its, origin, by = "OTU", all.x=TRUE)

its.30gen.names <- psz.its %>%
  filter(!is.na(Genus) & Genus != "NA" & !grepl("Incertae_sedis", Genus)) %>%
  group_by(Genus) %>%
  summarize(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 30) %>%
  pull(Genus) 
its.30gen.names
its.30gen <- subset(psz.its, psz.its$Genus %like% its.30gen.names & !is.na(SampleCategory))



#################################### MAKE PLOTS #################################### 

#get colors
#######
colors10 <- c(
  "#1f77b4", # Steel Blue
  "#ff7f0e", # Dark Orange
  "#2ca02c", # Forest Green
  "#9467bd", # Medium Purple
  "#e377c2", # Pink
  "#17becf", # Light Blue
  "#98df8a", # Light Green
  "#c5b0d5", # Light Purple
  "#9edae5", # Pale Blue
  "#393b79" # Dark Blue
)
colors20 <- c(
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

colors30 <- c(
  "#1f77b4", # Steel Blue
  "#ff7f0e", # Dark Orange
  "#2ca02c", # Forest Green
  "#9467bd", # Medium Purple
  "#d62728", # Brick Red
  "#8c564b", # Brown
  "#e377c2", # Pink
  "#7f7f7f", # Gray
  "#bcbd22", # Olive
  "#ff9896", # Light Salmon
  "#17becf", # Light Blue
  "#98df8a", # Light Green
  "#ffbb78", # Light Orange
  "#c5b0d5", # Light Purple
  "#c49c94", # Light Brown
  "#f7b6d2", # Light Pink
  "#dbdb8d", # Khaki
  "#9edae5", # Pale Blue
  "#393b79", # Dark Blue
  "#7b4173", # Dark Purple
  "#ffa500", # Orange
  "#ff9896", # Light Salmon
  "#00ced1", # Dark Turquoise
  "#32cd32", # Lime Green
  "#ba55d3", # Medium Orchid
  "#daa520", # Goldenrod
  "#00fa9a", # Medium Spring Green
  "#4682b4", # Steel Blue
  "#ff6347", # Tomato
  "#48d1cc"  # Medium Turquoise
)


######

its.soil.plot <- ggplot(its.soil.20gen, aes(y=Abundance, x=Tree, fill = Genus, color = Genus)) + 
  geom_bar(position="fill", stat="identity")+
  labs(x = "Tree species", y = "Relative abundance")+
  scale_fill_manual(values = colors20)+
  scale_color_manual(values = colors20)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13), axis.text=element_text(size=12))
its.soil.plot
ggsave("its_soil_stackedBarplot_top20genera.png", path = "Aus-Invasions-2023-Course/Aus23_ITS_Metabarcoding/StackedBarplots", width=8, height=10, dpi=300)

its.root.plot <- ggplot(its.root.20gen, aes(y=Abundance, x=Tree, fill = Genus, color = Genus)) + 
  geom_bar(position="fill", stat="identity")+
  labs(x = "Tree species", y = "Relative abundance")+
  scale_fill_manual(values = colors20)+ 
  scale_color_manual(values = colors20)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13), axis.text=element_text(size=12))
its.root.plot
ggsave("its_root_stackedBarplot_top20genera.png", path = "Aus-Invasions-2023-Course/Aus23_ITS_Metabarcoding/StackedBarplots", width=8, height=10, dpi=300)


cat.ord <- c("Eucalypt Root", "Pine Root", "Eucalypt Soil", "Pine Soil")
its.all.plot <- ggplot(its.30gen, aes(y=Abundance, x=factor(SampleCategory, cat.ord), fill = Genus, color = Genus)) + 
  geom_bar(position="fill", stat="identity")+
  labs(x = "Sample type", y = "Relative abundance")+
  scale_fill_manual(values = colors30)+ 
  scale_color_manual(values = colors30)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13), axis.text=element_text(size=12))
its.all.plot
ggsave("its_SoilRoot_stackedBarplot_top30genera.png", path = "Aus-Invasions-2023-Course/Aus23_ITS_Metabarcoding/StackedBarplots", width=9, height=8, dpi=300)

its.origin.plot <- ggplot(psz.its[!is.na(psz.its$SampleCategory),], aes(y=Abundance, x=factor(SampleCategory, cat.ord), fill = prevenance, color = prevenance)) + 
  geom_bar(position="fill", stat="identity")+
  labs(x = "Sample type", y = "Relative abundance")+
  scale_fill_manual(values = c("dodgerblue1", "chartreuse3", "gold"))+ 
  scale_color_manual(values = c("dodgerblue1", "chartreuse3", "gold"))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13), axis.text=element_text(size=12))
its.origin.plot
ggsave("its_native_introduced_stackedBarplot.png", path = "Aus-Invasions-2023-Course/Aus23_ITS_Metabarcoding/StackedBarplots", width=9, height=8, dpi=300)


psz.its.origin <- psz.its[!is.na(psz.its$SampleCategory) & !is.na(psz.its$prevenance) & !psz.its$prevenance == "uncertain",]
psz.its.origin <- psz.its.origin[order(psz.its.origin$prevenance),]
gen.ord <- unique(psz.its.origin$Genus)
gen.orig <- unique(psz.its.origin[c(9,12)])
origin.col <- c("#FAC898", "#FFA500", "#F28C28", "#E34A27", "#8B4000", "#4B0082", "#800080", "#8A2BE2", "#bc85fa", "#DA70D6", "#FF8BFF", "#DDA0DD", "#CEC2EB", "#E6E6FA")
origin.col2 <- c("#FAC898", "#FFA500", "#F28C28", "#E34A27", "#8B4000","#2e035e","#4B0082", "#5c07bc", "#8A2BE2", "#bc85fa", "#DA70D6", "#DDA0DD", "#CEC2EB", "#E6E6FA")
origin.col2 <- c("#FAC898", "#FFA500", "#F28C28", "#E34A27", "#8B4000","#2e035e","#4B0082", "#5c07bc", "#7409eb", "#8A2BE2", "#a75ef8", "#c08dfa", "#dabcfc", "#f4ebfe")


its.origin.only.plot <- ggplot(psz.its.origin, aes(y=Abundance, x=factor(SampleCategory, cat.ord), fill = factor(Genus, gen.ord), color = factor(Genus, gen.ord))) + 
  geom_bar(position="fill", stat="identity")+
  labs(x = "Sample type", y = "Relative abundance")+
  scale_fill_manual(values = origin.col2)+ 
  scale_color_manual(values = origin.col2)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13), axis.text=element_text(size=12))
its.origin.only.plot
ggsave("its_native_introduced_stackedBarplot.png", path = "Aus-Invasions-2023-Course/Aus23_ITS_Metabarcoding/StackedBarplots", width=9, height=8, dpi=300)



############################################################################################################
# EMF ONLY
############################################################################################################

# read in taxonomy list with functional guilds
its.tax.fxn <- read_tsv("Aus-Invasions-2023-Course/Aus23_ITS_Metabarcoding/ITS2_Dada2_repseqs97_taxonomy_edited_FunT_FunG.tsv")
colnames(its.tax.fxn)[1] <- "Feature_ID"
# get list of ectomycorrhizal OTUs
emf.otu.list <- its.tax.fxn %>%
  filter(grepl("Ectomycorrhizal", guild, ignore.case = TRUE)) %>%
  pull(Feature_ID)
# subset main df for only EMF taxa
psz.emf <- subset(psz.its, psz.its$OTU %in% emf.otu.list)

#Select top 10 most abundant genera
emf.10gen.names <- psz.emf %>%
  filter(!is.na(Genus) & Genus != "NA" & !grepl("Incertae_sedis", Genus)) %>%
  group_by(Genus) %>%
  summarize(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 10) %>%
  pull(Genus) 
emf.10gen.names
emf.10gen <- subset(psz.emf, psz.emf$Genus %like% emf.10gen.names & !is.na(SampleCategory))


cat.ord <- c("Eucalypt Root", "Pine Root", "Eucalypt Soil", "Pine Soil")
emf.plot <- ggplot(emf.10gen, aes(y=Abundance, x=factor(SampleCategory, cat.ord), fill = Genus, color = Genus)) + 
  geom_bar(position="fill", stat="identity")+
  labs(x = "Sample type", y = "Relative abundance")+
  scale_fill_manual(values = colors10)+ 
  scale_color_manual(values = colors10)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13), axis.text=element_text(size=12))
emf.plot
ggsave("EMF_SoilRoot_stackedBarplot_top10genera.png", path = "Aus-Invasions-2023-Course/Aus23_ITS_Metabarcoding/StackedBarplots", width=7, height=6, dpi=300)


tax.org <- subset(tax.test, tax.test$prevenance == NA)










