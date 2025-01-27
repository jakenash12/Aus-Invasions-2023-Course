#libraries
library(dplyr)

## Roots
# Read the CSV file into R (update file path if needed)
otu_data_AUS <- read.csv("~/Dropbox/NDSU_FOLDERS/Lennel/THESIS/NSF-IRES-AUSTRALIA/Aus-Invasions-2023-Course/Aus23_16S_Metabarcoding/OTUTables/OTU_16S_root_rel_melted.csv")

# Filter for Pine samples and get the top 25 OTUs based on abundance
top_25_pine_root <- otu_data_AUS %>%
  filter(Tree == "Pine") %>%   # Adjust the column and value to match your dataset
  group_by(OTU) %>%
  summarise(TotalAbundance = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 25)|>
  left_join(otu_data_AUS|> select(OTU,Kingdom:Species)|>unique())

# Filter for Euk samples and get the top 25 OTUs based on abundance
top_25_euk_root <- otu_data_AUS %>%
  filter(Tree == "Eucalypt") %>%    # Adjust the column and value to match your dataset
  group_by(OTU) %>%
  summarise(TotalAbundance = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 25)|>
  left_join(otu_data_AUS|> select(OTU,Kingdom:Species)|>unique())

## Soil

# Read the CSV file into R (update file path if needed)
otu_data_AUS_soil <- read.csv("~/Dropbox/NDSU_FOLDERS/Lennel/THESIS/NSF-IRES-AUSTRALIA/Aus-Invasions-2023-Course/Aus23_16S_Metabarcoding/OTUTables/OTU_16S_soil_rel_melted.csv")

# Filter for Pine samples and get the top 25 OTUs based on abundance
top_25_pine_soil <- otu_data_AUS_soil %>%
  filter(Tree == "Pine") %>%   # Adjust the column and value to match your dataset
  group_by(OTU) %>%
  summarise(TotalAbundance = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 25)|>
  left_join(otu_data_AUS_soil|> select(OTU,Kingdom:Species)|>unique())

# Filter for Euk samples and get the top 25 OTUs based on abundance
top_25_euk_soil <- otu_data_AUS_soil %>%
  filter(Tree == "Eucalypt") %>%    # Adjust the column and value to match your dataset
  group_by(OTU) %>%
  summarise(TotalAbundance = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 25)|>
  left_join(otu_data_AUS_soil|> select(OTU,Kingdom:Species)|>unique())

## ITS

## Roots
# Read the CSV file into R (update file path if needed)
otu_data_AUS_root_ITS <- read.csv("~/Dropbox/NDSU_FOLDERS/Lennel/THESIS/NSF-IRES-AUSTRALIA/Aus-Invasions-2023-Course/Aus23_ITS_Metabarcoding/OTUTables/OTU_ITS_root_rel_melted.csv")

# Filter for Pine samples and get the top 25 OTUs based on abundance
top_25_pine_root_ITS <- otu_data_AUS_root_ITS %>%
  filter(Tree == "Pine") %>%   # Adjust the column and value to match your dataset
  group_by(OTU) %>%
  summarise(TotalAbundance = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 25)|>
  left_join(otu_data_AUS_root_ITS|> select(OTU,Kingdom:Species)|>unique())

# Filter for Euk samples and get the top 25 OTUs based on abundance
top_25_euk_root_ITS <- otu_data_AUS_root_ITS %>%
  filter(Tree == "Eucalypt") %>%    # Adjust the column and value to match your dataset
  group_by(OTU) %>%
  summarise(TotalAbundance = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 25)|>
  left_join(otu_data_AUS_root_ITS|> select(OTU,Kingdom:Species)|>unique())

## Soil

# Read the CSV file into R (update file path if needed)
otu_data_AUS_soil_ITS <- read.csv("~/Dropbox/NDSU_FOLDERS/Lennel/THESIS/NSF-IRES-AUSTRALIA/Aus-Invasions-2023-Course/Aus23_ITS_Metabarcoding/OTUTables/OTU_ITS_soil_rel_melted.csv")

# Filter for Pine samples and get the top 25 OTUs based on abundance
top_25_pine_soil_ITS <- otu_data_AUS_soil_ITS %>%
  filter(Tree == "Pine") %>%   # Adjust the column and value to match your dataset
  group_by(OTU) %>%
  summarise(TotalAbundance = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 25)|>
  left_join(otu_data_AUS_soil_ITS|> select(OTU,Kingdom:Species)|>unique())

# Filter for Euk samples and get the top 25 OTUs based on abundance
top_25_euk_soil_ITS <- otu_data_AUS_soil_ITS %>%
  filter(Tree == "Eucalypt") %>%    # Adjust the column and value to match your dataset
  group_by(OTU) %>%
  summarise(TotalAbundance = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 25)|>
  left_join(otu_data_AUS_soil_ITS|> select(OTU,Kingdom:Species)|>unique())

## 18S 

## Roots
# Read the CSV file into R (update file path if needed)
otu_data_AUS_root_18S <- read.csv("~/Dropbox/NDSU_FOLDERS/Lennel/THESIS/NSF-IRES-AUSTRALIA/Aus-Invasions-2023-Course/Aus23_18S_Metabarcoding/OTUTables/OTU_18S_root_rel_melted.csv")

# Filter for Pine samples and get the top 25 OTUs based on abundance
top_25_pine_root_18S <- otu_data_AUS_root_18S %>%
  filter(Tree == "Pine") %>%   # Adjust the column and value to match your dataset
  group_by(OTU) %>%
  summarise(TotalAbundance = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 25)|>
  left_join(otu_data_AUS_root_18S|> select(OTU,Kingdom:Species)|>unique())

# Filter for Euk samples and get the top 25 OTUs based on abundance
top_25_euk_root_18s <- otu_data_AUS_root_18S %>%
  filter(Tree == "Eucalypt") %>%    # Adjust the column and value to match your dataset
  group_by(OTU) %>%
  summarise(TotalAbundance = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 25)|>
  left_join(otu_data_AUS_root_18S|> select(OTU,Kingdom:Species)|>unique())

## Soil

# Read the CSV file into R (update file path if needed)
otu_data_AUS_soil_18S <- read.csv("~/Dropbox/NDSU_FOLDERS/Lennel/THESIS/NSF-IRES-AUSTRALIA/Aus-Invasions-2023-Course/Aus23_18S_Metabarcoding/OTUTables/OTU_18S_soil_rel_melted.csv")

# Filter for Pine samples and get the top 25 OTUs based on abundance
top_25_pine_soil_18S <- otu_data_AUS_soil_18S %>%
  filter(Tree == "Pine") %>%   # Adjust the column and value to match your dataset
  group_by(OTU) %>%
  summarise(TotalAbundance = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 25)|>
  left_join(otu_data_AUS_soil_18S|> select(OTU,Kingdom:Species)|>unique())

# Filter for Euk samples and get the top 25 OTUs based on abundance
top_25_euk_soil_18S <- otu_data_AUS_soil_18S %>%
  filter(Tree == "Eucalypt") %>%   # Adjust the column and value to match your dataset
  group_by(OTU) %>%
  summarise(TotalAbundance = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 25)|>
  left_join(otu_data_AUS_soil_18S|> select(OTU,Kingdom:Species)|>unique())


# Create a list of data frames
list_of_tables_top_pine <- list(
  "16S.pine.soil" = top_25_pine_soil,
  "16S.pine.root" = top_25_pine_root,
  "ITS.pine.soil" = top_25_pine_soil_ITS,
  "ITS.pine.root" = top_25_pine_root_ITS,
  "18S.pine.soil" = top_25_pine_soil_18S,
  "18S.pine.root" = top_25_pine_root_18S
)

list_of_tables_top_euk <- list(
  "16S.euk.soil" = top_25_euk_soil,
  "16S.euk.root" = top_25_euk_root,
  "ITS.euk.soil" = top_25_euk_soil_ITS,
  "ITS.euk.root" = top_25_euk_root_ITS,
  "18S.euk.soil" = top_25_euk_soil_18S,
  "18S.euk.root" = top_25_euk_root_18s
)



# Write the list to an Excel file
writexl::write_xlsx(list_of_tables_top_pine, "~/Dropbox/NDSU_FOLDERS/Lennel/THESIS/NSF-IRES-AUSTRALIA/Aus-Invasions-2023-Course/Top25_OTUs_Pine.xlsx")
writexl::write_xlsx(list_of_tables_top_euk, "~/Dropbox/NDSU_FOLDERS/Lennel/THESIS/NSF-IRES-AUSTRALIA/Aus-Invasions-2023-Course/Top25_OTUs_Eucalypt.xlsx")


