# Make supplemental figure of litter % nitrogen by tree species
# 1/28/25 C. Vietorisz

library(tidyverse)

setwd("Aus-Invasions-2023-Course/")

# look at distribution of all litter C, N, P
mean(main$Pine_leafLitter_percN) # 0.56
mean(main$Euc_leafLitter_percN) # 0.69
hist(main$Pine_leafLitter_percN)
hist(main$Euc_leafLitter_percN)

mean(main$Pine_leafLitter_percC) # 50.8
mean(main$Euc_leafLitter_percC) # 53.3
hist(main$Pine_leafLitter_percC)
hist(main$Euc_leafLitter_percC)

mean(main$Pine_leafLitter_percP) # 0.024
mean(main$Euc_leafLitter_percP) # 0.015
hist(main$Pine_leafLitter_percP)
hist(main$Euc_leafLitter_percP)


# pivot dataframe longer
litN <- main[c("Pine_leafLitter_percN", "Euc_leafLitter_percN")]
colnames(litN) <- c("Pine", "Eucalypt")
litN_long <- pivot_longer(litN, cols = c("Pine", "Eucalypt"), names_to = "Litter_type", values_to = "LitterN")

# make plot
litNplot = ggplot(litN_long, aes(x=Litter_type, y=LitterN))+
  geom_point(size=1, colour= "darkgrey", position=position_jitter(width=0.1))+
  geom_boxplot(alpha = 0.35) +
  labs(x="Litter species", y="Litter % nitrogen")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=15),axis.title.y = element_text(size=15), axis.text=element_text(size=12), legend.key.size = unit(2, "lines"), legend.title = element_text(size=15), legend.text = element_text(size=12))
litNplot
ggsave("Soil_data/Euc_v_Pine_litter_percN.png", width = 3.3, height = 3.3, dpi = 300)

litN.lm <- lm(LitterN ~ Litter_type, data=litN_long)
anova(litN.lm)
plot(litN.lm)

