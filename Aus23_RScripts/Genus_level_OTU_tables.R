# Create genus-level OTU tables
# 1/6/24 CV

######## ITS ###########
# read in OTU tables with relative abundances
otu.its <- read.csv("Aus-Invasions-2023-Course/Aus23_ITS_Metabarcoding/OTUTables/OtuMatITS_rel.csv", row.names=1)

# read in taxonomy
its.tax <- read_tsv("Aus-Invasions-2023-Course/Aus23_ITS_Metabarcoding/ITS2_Dada2_repseqs97_taxonomy_edited_FunT_FunG.tsv")
its.tax.split <- its.tax %>% 
  separate(Kingdom, into = c("Phylum", "Class", "Order", "Family", "Genus"), sep = ";") %>%
  column_to_rownames(var = "Feature ID") %>%
  mutate(Kingdom = "Fungi") %>%
  select(Genus) %>%
  mutate(across(everything(), ~na_if(., "NA"))) %>%
  mutate(across(everything(), ~ifelse(grepl("Incertae_sedis", ., fixed = TRUE), NA, .)))

# adjust so OTU names are rows
otu.its.t <- as.data.frame(t(otu.its))
otu.its.t <- merge(otu.its.t, its.tax.split, by = 0, all.x=TRUE)
colnames(otu.its.t)[1] <- "OTU"
otu.its.t <- otu.its.t[!is.na(otu.its.t$Genus),]

# group by Genus
otu.genus.its <- otu.its.t %>%
  group_by(Genus) %>%
  summarise(across(where(is.numeric), sum)) %>%
  column_to_rownames("Genus") %>%
  t() %>%
  as.data.frame() 

# save as csv
write.csv(otu.genus.its, "Aus-Invasions-2023-Course/Aus23_ITS_Metabarcoding/OTUTables/OtuMatITS_rel_byGenus.csv")

######## Pooled data

# read in OTU tables with relative abundances
otu.its.pool <- read.csv("Aus-Invasions-2023-Course/Aus23_ITS_Metabarcoding/OTUTables/OtuMatITS_pooled_rel.csv", row.names=1) %>%
  filter(row.names(.) != "Neg_Control")

# adjust so OTU names are rows
otu.its.pool.t <- as.data.frame(t(otu.its.pool))
otu.its.pool.t <- merge(otu.its.pool.t, its.tax.split, by = 0, all.x=TRUE)
colnames(otu.its.pool.t)[1] <- "OTU"
otu.its.pool.t <- otu.its.pool.t[!is.na(otu.its.pool.t$Genus),]

# group by Genus
otu.genus.its.pool <- otu.its.pool.t %>%
  group_by(Genus) %>%
  summarise(across(where(is.numeric), sum)) %>%
  column_to_rownames("Genus") %>%
  t() %>%
  as.data.frame() 

# save as csv
write.csv(otu.genus.its.pool, "Aus-Invasions-2023-Course/Aus23_ITS_Metabarcoding/OTUTables/OtuMatITS_rel_byGenus_pooled.csv")



######## 16S ###########
# read in OTU tables with relative abundances
otu.16S <- read.csv("Aus-Invasions-2023-Course/Aus23_16S_Metabarcoding/OTUTables/OtuMat16S_rel.csv", row.names=1)

# read in taxonomy
bac.tax <- read_tsv("Aus-Invasions-2023-Course/Aus23_16S_Metabarcoding/Aus23_16S_SilvaTaxonomy_16S.tsv")[1:2]
bac.tax$Taxon <- str_replace_all(bac.tax$Taxon, "\\w__", "") # remove 
bac.tax.split <- bac.tax %>% 
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
  column_to_rownames(var = "Feature ID") %>%
  select(Genus) %>%
  mutate(across(everything(), ~na_if(., "NA"))) %>%
  mutate(across(everything(), ~ifelse(grepl("uncultured", ., fixed = TRUE), NA, .))) %>%
  mutate(across(everything(), ~ifelse(grepl("Unknown", ., fixed = TRUE), NA, .)))
bac.tax.split[] <- lapply(bac.tax.split, function(x) if(is.character(x)) trimws(x) else x) # remove the space at the beginning of each cell


# adjust so OTU names are rows
otu.16S.t <- as.data.frame(t(otu.16S))
otu.16S.t <- merge(otu.16S.t, bac.tax.split, by = 0, all.x=TRUE)
colnames(otu.16S.t)[1] <- "OTU"
otu.16S.t <- otu.16S.t[!is.na(otu.16S.t$Genus),]

# group by Genus
otu.genus.16S <- otu.16S.t %>%
  group_by(Genus) %>%
  summarise(across(where(is.numeric), sum)) %>%
  column_to_rownames("Genus") %>%
  t() %>%
  as.data.frame()

# save as csv
write.csv(otu.genus.16S, "Aus-Invasions-2023-Course/Aus23_16S_Metabarcoding/OTUTables/OtuMat16S_rel_byGenus.csv")

######## Pooled data

# read in OTU tables with relative abundances
otu.16S.pool <- read.csv("Aus-Invasions-2023-Course/Aus23_16S_Metabarcoding/OTUTables/OtuMat16S_pooled_rel.csv", row.names=1) %>%
  filter(row.names(.) != "Neg_Control")

# adjust so OTU names are rows
otu.16S.pool.t <- as.data.frame(t(otu.16S.pool))
otu.16S.pool.t <- merge(otu.16S.pool.t, bac.tax.split, by = 0, all.x=TRUE)
colnames(otu.16S.pool.t)[1] <- "OTU"
otu.16S.pool.t <- otu.16S.pool.t[!is.na(otu.16S.pool.t$Genus),]

# group by Genus
otu.genus.16S.pool <- otu.16S.pool.t %>%
  group_by(Genus) %>%
  summarise(across(where(is.numeric), sum)) %>%
  column_to_rownames("Genus") %>%
  t() %>%
  as.data.frame() 

# save as csv
write.csv(otu.genus.16S.pool, "Aus-Invasions-2023-Course/Aus23_16S_Metabarcoding/OTUTables/OtuMat16S_rel_byGenus_pooled.csv")



