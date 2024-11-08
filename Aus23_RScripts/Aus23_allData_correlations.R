# Univariate regressions between microbial metrics and soil nutrient, biolog, and litter nutrient data
# 7 Nov 2024
# Corinne V

setwd("~/Desktop/GitHub/Aus-Invasions-2023-Course")

# read in data
main <- read.csv("Merged_data/Aus23_allData_7Nov24.csv", row.names=1) # pooled
CNunpool <- read.csv("Merged_data/Aus23_microbes_soilCN_unpooled.csv", row.names=1) # pooled

########### make a correlation heatmap between all variables ###############
library(Hmisc)
#select out continuous variables only
corr_vars <- main[c(4:7,9:77,80:205)]
r <- rcorr(as.matrix(corr_vars))$r
p <- rcorr(as.matrix(corr_vars))$P
#save correlations and p values
write.csv(r, "Correlations/All_variable_correlations.csv")
write.csv(p, "Correlations/All_variable_correlation_pvals.csv")

