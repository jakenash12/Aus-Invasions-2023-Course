#BB
#10/30/24
#Trying some regressions to begin with

library(ggplot2)

ds <- read.csv("Aus23_CNP_pooled_biolog.csv")


# Gadgil effect -----------------------------------------------------------


mod <- lm(N.Acetyl.D.Glucosamine ~ ECM_richness_soil + perc_N + ergosterol + TreeSpecies, ds)
summary(mod)

mod <- lm(ergosterol ~ ECM_abund_soil + perc_N + ergosterol + TreeSpecies, ds)
summary(mod)
#interesting...

#link who they are and what they do

#manylm. Try this function