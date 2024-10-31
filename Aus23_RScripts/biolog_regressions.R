#BB
#10/30/24
#Trying some regressions to begin with

library(ggplot2)
library(mvabund)

ds <- read.csv("Aus23_CNP_pooled_biolog.csv")


# Gadgil effect -----------------------------------------------------------
#How does ECM affect litter decomp in the two types of soils?

mod <- lm(N.Acetyl.D.Glucosamine ~ ECM_richness_soil + perc_N + ergosterol + TreeSpecies, ds)
summary(mod)

mod <- lm(ergosterol ~ ECM_abund_soil + perc_N + ergosterol + TreeSpecies, ds)
summary(mod)
#interesting... mpisture should be in this too

# many <- manylm(N.Acetyl.D.Glucosamine~.,data=ds)
# many



#link who they are and what they do

#manylm. Try this function

#Do some sort of AIC-based model selection with all predictors
#Response variables; Total CNP, decomp, 
#Decomp 