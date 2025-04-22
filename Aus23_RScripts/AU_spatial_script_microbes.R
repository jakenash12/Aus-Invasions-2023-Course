#AU Manuscript Analyses: Spatial Analysis
#Microbial; in response to manuscript reviewer's comment

#Libraries
library(spdep)
library(sf)
library(dplyr)
library(ggplot2)
library(vegan)
library(spgwr)
library(gstat)
library(viridis)

#Working Directory
setwd("/Users/crvietorisz/Desktop/BU/PhD/Pine invasion/Australia NSF course/Belanglo project/Rscripts/autocorrelation")

#Dataset
coord <- read.csv("AUlatlong.csv")
all_data <- read.csv("/Users/crvietorisz/Desktop/GitHub/Aus-Invasions-2023-Course/Merged_data/Aus23_allData_18Jan25.csv")

sigvar.order <-c("treeName", "Tween_40_BiologDay5", "L.Serine_BiologDay5", 
                 "Glycyl.L.Glutamic_Acid_BiologDay5",
                 "i.Erythritol_BiologDay5", "Putrescine_BiologDay5", "L.Asparagine_BiologDay5", "soil_moisture", "ergosterol", 
                 "Litter_OLayer_ergosterol", "Litter_avg_ergosterol",
                 "pine_litter_prop", "euc_litter_prop", "Euc_leafLitter_percP", "Bac_Shannon_soil")
sub_data <- all_data[sigvar.order]


#Ensure latitude and longitude values are numeric
data <- data %>%
  mutate(
    latitude = as.numeric(latitude),
    longitude = as.numeric(longitude)
  )

#Convert the dataset to an sf object
data_sf <- st_as_sf(data, coords = c("longitude", "latitude"), crs = 4326)

#Create a spatial weights matrix
coords <- st_coordinates(data_sf)
knn <- knearneigh(coords, k = 4)  #Examining 4 nearest neighbors
knn_weights <- nb2listw(knn2nb(knn), style = "W")

#Define groups of variables
#ITS Soil
#variables <- c("Phialocephala", "Thelephora", "Athelopsis", "Diplodia", "Ascocorticium", "Tomentella", "Leohumicola", "Clavulinopsis", "Pseudocamarosporium", "Flavocillium")

#ITS Roots
#variables <- c("Phialocephala", "Talaromyces", "Devriesia", "Aspergillus", "Penicillium", "Auricularia", "Pezoloma", "Cryptosporiopsis", "Pezicula", "Mollisia")

#16S Soil
#variables <- c("Gryllotalpicola", "Limnobacter", "Acidisoma", "Dinghuibacter", "Acidiphilium", "11-24", "SM1A02", "Subgroup_5", "A4b", "Mesorhizobium")

#16S Roots
variables <- c("Lechevalieria", "Nocardia", "Acidicaldus", "Nevskia", "Occallatibacter", "SBR1031", "Novosphingobium", "Collimonas", "A4b", "Bryocella")

######--- 1. Moran's I for Global Spatial Autocorrelation ---######
results_global <- list()
for (var in variables) {
  if (var %in% colnames(data)) {
    moran <- moran.test(data[[var]], listw = knn_weights)
    results_global[[var]] <- list(
      variable = var,
      moran_I = moran$estimate[1],
      p_value = moran$p.value
    )
  }
}

#Convert results to a data frame
results_global_df <- do.call(rbind, lapply(results_global, as.data.frame))
results_global_df <- as.data.frame(results_global_df)

#Save results
write.csv(results_global_df, "OtuMat16S_rare_Roots_Top5Genus_pooled_GlobalMorans.csv", row.names = FALSE)
print(results_global_df)

