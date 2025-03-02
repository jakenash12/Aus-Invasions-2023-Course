library(tidyverse)
library(magrittr)
library(readxl)

#reads in curated sheet of native and introduced fungi
#also creates a new Species column that gives duplicates
#(i.e. multiple Laccara or Tomentella sp) unique numbers
ITS2_Curated_NativeIntroduced=
  read_xlsx("../Aus23_ITS_Metabarcoding/ITS2_Curated_NativeIntroduced.xlsx") %>%
  group_by(Species) %>%
  arrange(Species, desc(Provenance == "Introduced"), .by_group = TRUE) %>%
  mutate(Species_unique=case_when(Species %in% c("Laccaria_sp", "Tomentella_sp") ~ paste(Species, row_number(), sep = " "),
                                  .default=Species)) %>%
  ungroup()

#generates more curated list of fungi of interest
SelectFungi=c("Suillus_quiescens", "Suillus_luteus",
              "Rhizopogon_pseudoroseolus", "Rhizopogon_verii",
              "Rhizopogon_evadens", "Lactarius_deliciosus",
              "Amanita_muscaria", "Thelephora_terrestris",
              "Phialocephela_sp", 
              )

OtuMatITS_pooled_rel_curated=
  OtuMatITS_pooled_rel %>%
  select(ITS2_Curated_NativeIntroduced$`Feature ID`) %>%
  mutate(TreeSampleType=rownames(.)) %>% 
  gather("otu","Abundance",-TreeSampleType) %>%
  left_join(Metadata_ITS_pooled) %>%
  filter(SampleType!="Neg") %>%
  left_join(ITS2_Curated_NativeIntroduced, by=c("otu"="Feature ID")) %>%
  mutate(Species = gsub("_", " ", Species)) %>%
  mutate(Species = gsub(" sp", " sp.", Species)) %>%
  mutate(TreeSpecies = gsub("Eucalyptus", "Eucalypt", TreeSpecies))

Introduced_plot =
  OtuMatITS_pooled_rel_curated %>%
  filter(Provenance == "Introduced") %>%
  mutate(otu = factor(otu, levels = unique(c(
    otu[Genus == "Suillus"], 
    otu[Genus == "Rhizopogon"], 
    otu[Genus == "Phialocephala"],
    otu[!Genus %in% c("Suillus", "Rhizopogon", "Phialocephala")]
  )))) %>%
  ggplot(aes(x = SampleType,
             y = Abundance,
             color = TreeSpecies)) +
  geom_boxplot(
    position = position_dodge(width = 0.8), # Increased width for more spacing
    outlier.shape = NA
  ) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.2, # Jittering amount on x-axis
      dodge.width = 0.8   # Dodge width to align with boxplots
    ),
    size = 2.5,
    alpha = 0.6
  ) +
  facet_wrap(
    . ~ otu, 
    scales = "free_y", 
    labeller = labeller(otu = function(otu_levels) {
      OtuMatITS_pooled_rel_curated$Species_unique[match(otu_levels, OtuMatITS_pooled_rel_curated$otu)]
    }),
    strip.position = "top"
  ) +
  scale_color_manual(values = c("#FF9900", "#000DCC")) +
  theme_test() +
  theme(
    axis.text = element_text(colour = "black"),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(),
    strip.text = element_text(margin = margin(b = 5, t = 5)) # Adds padding to facet strip
  ) +
  ggtitle("Introduced Fungi") +
  theme(plot.title = element_text(hjust = 0.5))

Native_plot =
  ggplot(
    filter(OtuMatITS_pooled_rel_curated, Provenance == "Native"), 
    aes(x = SampleType,
        y = Abundance,
        color = TreeSpecies)
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.8), # Increased width for more spacing
    outlier.shape = NA
  ) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.2, # Jittering amount on x-axis
      dodge.width = 0.8   # Dodge width to align with boxplots
    ),
    size = 2.5,
    alpha = 0.6
  ) +
  facet_wrap(
    . ~ otu, 
    scales = "free_y", 
    labeller = labeller(otu = function(otu_levels) {
      OtuMatITS_pooled_rel_curated$Species_unique[match(otu_levels, OtuMatITS_pooled_rel_curated$otu)]
    }),
    strip.position = "top"
  ) +
  scale_color_manual(values = c("#FF9900", "#000DCC")) +
  theme_test() +
  theme(
    axis.text = element_text(colour = "black"),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(),
    strip.text = element_text(margin = margin(b = 5, t = 5)) # Adds padding to facet strip
  ) +
  ggtitle("Native Australian Fungi") +
  theme(plot.title = element_text(hjust = 0.5))


