#### Funal Traits and Fun Guild combo
### Author: Lennel A. Camuy-Velez

## FungalTraits
#install.packages("readxl")

## Fungal Traits DB

FT.db <- readxl::read_excel("~/Downloads/13225_2020_466_MOESM4_ESM.xlsx")

## Cleaning the taxonomic ranks:

clean_taxonomic_rank <- function(rank_column) {
  sub("^.*__", "", rank_column)
}

# Define a function to fill traits based on hierarchical taxonomy
fill_traits <- function(row, trait_data) {
  traits <- list(primary_lifestyle = NA, Secondary_lifestyle = NA)
  
  if (!is.na(row$Genus)) {
    matched_row <- trait_data %>% filter(GENUS == row$Genus)
    if (nrow(matched_row) > 0) {
      traits$primary_lifestyle <- matched_row$primary_lifestyle[1]
      traits$Secondary_lifestyle <- matched_row$Secondary_lifestyle[1]
      return(traits)
    }
  }
  
  if (is.na(traits$primary_lifestyle) && !is.na(row$Family)) {
    matched_row <- trait_data %>% filter(Family== row$Family)
    if (nrow(matched_row) > 0) {
      traits$primary_lifestyle <- matched_row$primary_lifestyle[1]
      traits$Secondary_lifestyle <- matched_row$Secondary_lifestyle[1]
      return(traits)
    }
  }
  
  if (is.na(traits$primary_lifestyle) && !is.na(row$Order)) {
    matched_row <- trait_data %>% filter(Order == row$Order)
    if (nrow(matched_row) > 0) {
      traits$primary_lifestyle <- matched_row$primary_lifestyle[1]
      traits$Secondary_lifestyle <- matched_row$Secondary_lifestyle[1]
      return(traits)
    }
  }
  
  if (is.na(traits$primary_lifestyle) && !is.na(row$Class)) {
    matched_row <- trait_data %>% filter(Class == row$Class)
    if (nrow(matched_row) > 0) {
      traits$primary_lifestyle <- matched_row$primary_lifestyle[1]
      traits$Secondary_lifestyle <- matched_row$Secondary_lifestyle[1]
      return(traits)
    }
  }
  if (is.na(traits$primary_lifestyle) && !is.na(row$Phylum)) {
    matched_row <- trait_data %>% filter(Phylum == row$Phylum)
    if (nrow(matched_row) > 0) {
      traits$primary_lifestyle <- matched_row$primary_lifestyle[1]
      traits$Secondary_lifestyle <- matched_row$Secondary_lifestyle[1]
      return(traits)
    }
  }
  
  return(traits)
}

Australia<-readxl::read_excel("~/Dropbox/NDSU_FOLDERS/Lennel/THESIS/NSF-IRES-AUSTRALIA/ITS2_Dada2_repseqs97_taxonomy_edited.xlsx")


# Apply the function to fill primary_lifestyle and Secondary_lifestyle
filled_traits.Australia <- Australia %>%
  rowwise() %>%
  mutate(traits = list(fill_traits(cur_data(), FT.db))) %>%
  unnest_wider(traits) %>%
  ungroup()

## FunGuild
library(FUNGuildR)
FUNGuild.db <- get_funguild_db()


filled_traits.Australia |>
  unite(col = "Kingdom","Phylum","Class","Order","Family","Genus",sep = ";")|>
  funguild_assign(db = FUNGuild.db,tax_col = "Kingdom") -> filled_traits.guilds.Australia

