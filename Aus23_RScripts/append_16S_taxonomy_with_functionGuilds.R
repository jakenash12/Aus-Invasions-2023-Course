# Append 16S Silva taxonomy table with functional guilds
# Using functional guild database created by Corinne's labmate Zoey: https://github.com/zoey-rw/soil_bacteria_functional_groups/tree/main
# script by C. Vietorisz, last edited 9/27/24 CV

library(tidyverse)

silva.tax <- read.delim("/Users/moniquegagnon/Desktop/BU/PhD/Pine invasion/Australia NSF course/Belanglo project/Data/Soil DNA/Taxonomy/Aus23_16S_SilvaTaxonomy_16S.tsv",sep="\t")

#import taxonomy table and parse into taxonomic levels
extract_text <- function(x) {
  sub(".+__", "", x)
}
silva.tax=
  read_tsv("/Users/moniquegagnon/Desktop/BU/PhD/Pine invasion/Australia NSF course/Belanglo project/Data/Soil DNA/Taxonomy/Aus23_16S_SilvaTaxonomy_16S.tsv") %>%
  rename("otu"="Feature ID") %>%
  select(-Confidence) %>%
  separate_wider_delim(too_few="align_start", cols = Taxon, delim = ";", names = c("Kingdom", "Phylum","Class", "Order","Family","Genus","Species")) %>%
  mutate(across(!contains("otu"), ~extract_text(.)))

#read in database
bac_fun <- read.csv("/Users/moniquegagnon/Desktop/BU/PhD/Databases/Bac_fxn_groups_database_Werbin.csv") #functional group .csv from Katie 

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

#combine with your taxonomy table
final_asv_fxn <- merge(silva.tax, fun_genus, by = "Genus", all.x = TRUE)
final_asv_fxn <- merge(final_asv_fxn, fun_family, by = "Family", all.x = TRUE)
final_asv_fxn <- merge(final_asv_fxn, fun_order, by = "Order", all.x = TRUE)
final_asv_fxn <- merge(final_asv_fxn, fun_class, by = "Class", all.x = TRUE)
final_asv_fxn <- merge(final_asv_fxn, fun_phylum, by = "Phylum", all.x = TRUE)

# create a new column that combines all possible functional categories 
final_asv_fxn$Functional_group <- paste(final_asv_fxn$system_g, final_asv_fxn$system_f, final_asv_fxn$system_o, final_asv_fxn$system_c, final_asv_fxn$system_p)
final_asv_fxn$Functional_group <- gsub("NA ","",final_asv_fxn$Functional_group)
final_asv_fxn$Functional_group <- gsub(" NA","",final_asv_fxn$Function)
# create a new column that combines all possible functions 
final_asv_fxn$Function <- paste(final_asv_fxn$function_g, final_asv_fxn$function_f, final_asv_fxn$function_o, final_asv_fxn$function_c, final_asv_fxn$function_p)
final_asv_fxn$Function <- gsub("NA ","",final_asv_fxn$Function)
final_asv_fxn$Function <- gsub(" NA","",final_asv_fxn$Function)

#Create a column that identifies only if the taxa is a copiotroph, oligotroph, or neither
final_asv_fxn <- final_asv_fxn %>%
  mutate(Copio_oligo = case_when(
    grepl("Copiotroph", Function) ~ "Copiotroph",
    grepl("Oligotroph", Function) ~ "Oligotroph",
    TRUE ~ NA_character_
  ))

write.csv(final_asv_fxn, file = "/Users/moniquegagnon/Desktop/BU/PhD/Pine invasion/Australia NSF course/Belanglo project/Data/Soil DNA/Taxonomy/Aus23_16S_SilvaTaxonomy_FunctionalGuilds.csv")

