library(tidyverse)
library(magrittr)
library(readxl)

#reads in curated sheet of native and introduced fungi
ITS2_Curated_NativeIntroduced=
  read_xlsx("../Aus23_ITS_Metabarcoding/ITS2_Curated_NativeIntroduced.xlsx")

OtuMatITS_pooled_rel_curated=
  OtuMatITS_pooled_rel %>%
  select(ITS2_Curated_NativeIntroduced$`Feature ID`) %>%
  mutate(TreeSampleType=rownames(.)) %>%
  gather("otu","Abundance",-TreeSampleType) %>%
  left_join(Metadata_ITS_pooled) %>%
  filter(SampleType!="Neg") %>%
  left_join(ITS2_Curated_NativeIntroduced, by=c("otu"="Feature ID"))

Introduced_plot=
  OtuMatITS_pooled_rel_curated %>%
  filter(Provenance=="Introduced") %>%
  mutate(otu = factor(otu, levels = unique(c(
    otu[Genus == "Suillus"], 
    otu[Genus == "Rhizopogon"], 
    otu[!Genus %in% c("Suillus", "Rhizopogon")]
  )))) %>%
  ggplot(aes(x=SampleType,
             y=Abundance,
             color=TreeSpecies)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    outlier.shape = NA
  ) + # Align boxplots properly for grouped categories
  geom_jitter(
    position = position_dodge(width = 0.75),
    size = 3,
    alpha = 0.6
  ) +
  facet_wrap(.~otu, scales = "free_y", labeller = labeller(otu = function(otu_levels) {
    OtuMatITS_pooled_rel_curated$Species[match(otu_levels, OtuMatITS_pooled_rel_curated$otu)]
  })) +
  scale_color_manual(values = c("#FF9900", "#000DCC")) +
  theme(axis.text = element_text(colour="black"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  ggtitle("Introduced Fungi") +
  theme(plot.title = element_text(hjust = 0.5))

Native_plot=
  ggplot(filter(OtuMatITS_pooled_rel_curated, Provenance=="Native"), 
         aes(x=SampleType,
             y=Abundance,
             color=TreeSpecies)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    outlier.shape = NA
  ) + # Align boxplots properly for grouped categories
  geom_jitter(
    position = position_dodge(width = 0.75),
    size = 3,
    alpha = 0.6
  ) +
  facet_wrap(.~otu, scales = "free_y", labeller = labeller(otu = function(otu_levels) {
    OtuMatITS_pooled_rel_curated$Species[match(otu_levels, OtuMatITS_pooled_rel_curated$otu)]
  })) +
  scale_color_manual(values = c("#FF9900", "#000DCC")) +
  theme(axis.text = element_text(colour="black"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  ggtitle("Native Australian Fungi") +
  theme(plot.title = element_text(hjust = 0.5))

