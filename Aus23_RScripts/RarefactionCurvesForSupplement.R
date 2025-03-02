#Generates rarefaction curves that are included in supplemental file
#color codes curves by tree species and sample type

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")

RareCurve_16S_plot=
  RareCurve_16S %>%
  rename(SampleID = Site) %>%
  left_join(Metadata_16S)

Rarecurve_ggplot_16S=
  ggplot(filter(RareCurve_16S_plot, SampleType!="Neg"), aes(x=Sample, y=Species, group=SampleID, color=TreeSpecies_SampleType)) +
  geom_line(size=1) +
  geom_vline(xintercept=rare_depth_16S, linetype=2) +
  ylab("# ASVs") +
  xlab("# Sequences") +
  scale_colour_manual(values=cbPalette) +
  theme_test() +
  ggtitle("16S Rarefaction curve - Depth=24860") +
  theme(legend.position="none")

RareCurve_ITS_plot=
  RareCurve_ITs %>%
  rename(SampleID = Site) %>%
  left_join(Metadata_16S)

Rarecurve_ggplot_ITS=
  ggplot(filter(RareCurve_ITS_plot, SampleType!="Neg"), aes(x=Sample, y=Species, group=SampleID, color=TreeSpecies_SampleType)) +
  geom_line(size=1) +
  geom_vline(xintercept=rare_depth_ITS, linetype=2) +
  ylab("# OTUs") +
  xlab("# Sequences") +
  scale_colour_manual(values=cbPalette) +
  theme_test() +
  ggtitle("ITS Rarefaction curve - Depth=70218") +
  theme(legend.position="none")

RareCurve_18S_plot=
  RareCurve_18S %>%
  rename(SampleID = Site) %>%
  left_join(Metadata_16S)

Rarecurve_ggplot_18S=
  ggplot(filter(RareCurve_18S_plot, SampleType!="Neg"), aes(x=Sample, y=Species, group=SampleID, color=TreeSpecies_SampleType)) +
  geom_line(size=1) +
  geom_vline(xintercept=rare_depth_18S, linetype=2) +
  ylab("# ASVs") +
  xlab("# Sequences") +
  scale_colour_manual(values=cbPalette) +
  theme_test() +
  ggtitle("18S Rarefaction curve - Depth=2053") +
  theme(legend.position="none")

plot_grid(Rarecurve_ggplot_ITS,
          Rarecurve_ggplot_18S,
          Rarecurve_ggplot_16S, ncol=3)

 View(RareCurve_16S)
