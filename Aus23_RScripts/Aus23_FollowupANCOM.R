ANCOM_ITS_pooled_soil_df %>%
  dplyr::filter(Genus=="Rhizopogon") %>%
  View
Phyloseq_pooled_ITS_soil %>%
  otu_table %>%
  as.data.frame %>% 
  mutate(TreeSampleType=rownames(.)) %>%
  left_join(Metadata_ITS_pooled) %>%
  ggplot(aes(x=TreeSpecies, y=f90e77a4f856e45f07957b2b61ef8ef3)) +
  geom_boxplot()
