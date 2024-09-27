# White pine 2021: soil DNA 16S
# Add functional annotations to bacterial data using Zoey Werbin's literature database: https://github.com/zoey-rw/soil_bacteria_functional_groups/tree/main
# 7/10/23

#Add functional groups

#read in 16S ASV table
asv <- read.csv("ASV.csv") #change to your ASV table

#read in database
bac_fun_orig <- read.csv("Bac_fxn_groups_database_Werbin.csv") #functional group .csv from Katie 
bac_fun <- read.csv("Bac_fxn_groups_database_Werbin.csv") #functional group .csv from Katie 

#this is super clunky im sorry but this picks out the taxonomic level that has the functional assignment, then merges that with the full df based on that taxonomic level
fun_genus <- bac_fun[bac_fun$Taxonomic.level == "Genus",]
fun_genus <- fun_genus[2:4]
colnames(fun_genus) <- c("Genus", "system_g", "function_g")

fun_family <- bac_fun[bac_fun$Taxonomic.level == "Family",]
fun_family <- fun_family[2:4]
colnames(fun_family) <- c("Family", "system_f", "function_f")

fun_order <- bac_fun[bac_fun$Taxonomic.level == "Order",]
fun_order <- fun_order[2:4]
colnames(fun_order) <- c("Order", "system_o", "function_o")

fun_class <- bac_fun[bac_fun$Taxonomic.level == "Class",]
fun_class <- fun_class[2:4]
colnames(fun_class) <- c("Class", "system_c", "function_c")

fun_phylum <- bac_fun[bac_fun$Taxonomic.level == "Phylum",]
fun_phylum <- fun_phylum[2:4]
colnames(fun_phylum) <- c("Phylum", "system_p", "function_p")

final_asv_fxn <- merge(asv, fun_genus, by = "Genus", all.x = TRUE)
final_asv_fxn <- merge(final_asv_fxn, fun_family, by = "Family", all.x = TRUE)
final_asv_fxn <- merge(final_asv_fxn, fun_order, by = "Order", all.x = TRUE)
final_asv_fxn <- merge(final_asv_fxn, fun_class, by = "Class", all.x = TRUE)
final_asv_fxn <- merge(final_asv_fxn, fun_phylum, by = "Phylum", all.x = TRUE)

final_asv_fxn$Functional_group <- paste(final_asv_fxn$system_g, final_asv_fxn$system_f, final_asv_fxn$system_o, final_asv_fxn$system_c, final_asv_fxn$system_p)
final_asv_fxn$Functional_group <- gsub("NA ","",final_asv_fxn$Functional_group)
final_asv_fxn$Functional_group <- gsub(" NA","",final_asv_fxn$Function)

final_asv_fxn$Function <- paste(final_asv_fxn$function_g, final_asv_fxn$function_f, final_asv_fxn$function_o, final_asv_fxn$function_c, final_asv_fxn$function_p)
final_asv_fxn$Function <- gsub("NA ","",final_asv_fxn$Function)
final_asv_fxn$Function <- gsub(" NA","",final_asv_fxn$Function)

write.csv(final_asv_fxn, file = "/Users/moniquegagnon/Desktop/BU/PhD/White pine/Soil DNA/ASV_tables/16S/ASV_functions_16S_rationorm_psz.csv")


