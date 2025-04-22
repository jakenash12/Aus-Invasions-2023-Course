#examining microbial species composition using PCoA axis 1 scores
library(car)
library(emmeans)
pco.scores.autocorr=read.csv("/Users/crvietorisz/Desktop/BU/PhD/Pine invasion/Australia NSF course/Belanglo project/Rscripts/autocorrelation/pco.scores.autocorr.csv", header = T)
#Does microbial community composition vary with location?
lm.result.fung.soil.pcoa1=lm(fung.soil.pcoa1 ~ latitude	+ longitude, data = pco.scores.autocorr)
Anova(lm.result.fung.soil.pcoa1)
summary(lm.result.fung.soil.pcoa1)
#Fungal composition in the soil does not vary significantly with location.
lm.result.bact.soil.pcoa1=lm(bact.soil.pcoa1 ~ latitude	+ longitude, pco.scores.autocorr)
Anova(lm.result.bact.soil.pcoa1)
summary(lm.result.bact.soil.pcoa1)
#Bacterial composition in the soil does not vary significantly with location.
lm.result.fung.root.pcoa1=lm(fung.root.pcoa1 ~ latitude	+ longitude, pco.scores.autocorr)
Anova(lm.result.fung.root.pcoa1)
summary(lm.result.fung.root.pcoa1)
#Fungal composition associated with roots does not vary significantly with location.
lm.result.bact.root.pcoa1=lm(bact.root.pcoa1 ~ latitude	+ longitude, pco.scores.autocorr)
Anova(lm.result.bact.root.pcoa1)
summary(lm.result.bact.root.pcoa1)
#Bacterial composition associated with roots does not vary significantly with location.
lm.result.amf.root.pcoa1=lm(amf.root.pcoa1 ~ latitude	+ longitude, pco.scores.autocorr)
Anova(lm.result.amf.root.pcoa1)
summary(lm.result.amf.root.pcoa1)
#AMF composition associated with roots does not vary significantly with location.
lm.result.amf.soil.pcoa1=lm(amf.soil.pcoa1 ~ latitude	+ longitude, pco.scores.autocorr)
Anova(lm.result.amf.soil.pcoa1)
summary(lm.result.amf.soil.pcoa1)
#AMF composition in the soil does not vary significantly with location.

#How does microbial community composition vary with respect to tree species and dbh after accounting for location using lat-long coordinates as covariates?
lm.result.fung.soil.pcoa1=lm(fung.soil.pcoa1 ~ latitude	+ longitude+DBH.cm*Focal_tree_spp, data = pco.scores.autocorr)
Anova(lm.result.fung.soil.pcoa1)
species.means=emmeans(lm.result.fung.soil.pcoa1, ~Focal_tree_spp)
species.means
summary(lm.result.fung.soil.pcoa1)
#Fungal composition in the soil does not vary significantly in response to tree species, dbh, or the interaction.
lm.result.bact.soil.pcoa1=lm(bact.soil.pcoa1 ~ latitude	+ longitude+DBH.cm*Focal_tree_spp, pco.scores.autocorr)
Anova(lm.result.bact.soil.pcoa1)
species.means=emmeans(lm.result.bact.soil.pcoa1, ~Focal_tree_spp)
species.means
summary(lm.result.bact.soil.pcoa1)
#Bacterial composition in the soil does not vary significantly in response to tree species, dbh, or the interaction.
lm.result.fung.root.pcoa1=lm(fung.root.pcoa1 ~ latitude	+ longitude+DBH.cm*Focal_tree_spp, pco.scores.autocorr)
Anova(lm.result.fung.root.pcoa1)
species.means=emmeans(lm.result.fung.root.pcoa1, ~Focal_tree_spp)
species.means
summary(lm.result.fung.root.pcoa1)
#Fungal composition associated with roots varies significantly in response to tree species and dbh, but not the interaction. Pines are associated with positive scores; Eucs with negative scores. Scores become increasingly positive with increasing dbh.
lm.result.bact.root.pcoa1=lm(bact.root.pcoa1 ~ latitude	+ longitude+DBH.cm*Focal_tree_spp, pco.scores.autocorr)
Anova(lm.result.bact.root.pcoa1)
species.means=emmeans(lm.result.bact.root.pcoa1, ~Focal_tree_spp)
species.means
summary(lm.result.bact.root.pcoa1)
#Bacterial composition associated with roots varies significantly in response to dbh. Scores become increasingly positive with increasing dbh.
lm.result.amf.root.pcoa1=lm(amf.root.pcoa1 ~ latitude	+ longitude+DBH.cm*Focal_tree_spp, pco.scores.autocorr)
Anova(lm.result.amf.root.pcoa1)
species.means=emmeans(lm.result.amf.root.pcoa1, ~Focal_tree_spp)
species.means
summary(lm.result.amf.root.pcoa1)
#AMF composition associated with roots does not vary significantly in response to tree species, dbh, or the interaction.
lm.result.amf.soil.pcoa1=lm(amf.soil.pcoa1 ~ latitude	+ longitude+DBH.cm*Focal_tree_spp, pco.scores.autocorr)
Anova(lm.result.amf.soil.pcoa1)
species.means=emmeans(lm.result.amf.soil.pcoa1, ~Focal_tree_spp)
species.means
summary(lm.result.amf.soil.pcoa1)
#AMF composition associated with soil does vary significantly in response to dbh, but not species. Scores becoming increasingly negative with increasing DBH.
library(ade4)
# make an object that is a list of the root-associated fungal PCoA scores.
root.fung.scores=pco.scores.autocorr[-c(1:5,7:23)]
head(root.fung.scores)
# make an object that is a list of the lat-long coordinates.
lat.long.coord=pco.scores.autocorr[-c(1:21)]
head(lat.long.coord)
#make a list of distances among the root.fung.scores
root.fung.dist=dist(root.fung.scores)
#convert to a matrix. Note that it is not necessary to do this step. It's for visualization.
root.fung.mat=as.matrix(root.fung.dist)
head(root.fung.mat)
#make a list of distances from the lat-long coordinates
lat.long.dist=dist(lat.long.coord)
#convert to a matrix
lat.long.mat=as.matrix(lat.long.dist)
head(lat.long.mat)
root.fung.test = mantel.rtest(lat.long.dist, root.fung.dist, nrepet = 9999)
root.fung.test
