library(tidyverse)
library(googledrive)
library(stringr)
library(magrittr)
library(ggplot2)
library(rmarkdown)

#Reads in substrate plate map - fill in path to where you downloaded file
sub_temp_path=''
substrate_template=
  read.csv(paste(sub_temp_path, 'Ecoplate_substrate_template.csv', sep=''), header=FALSE) %>%
  as.matrix

#creates a vector with all the subsrate names that will be looped through
substrate_vector=
  substrate_template %>%
  as.character %>% 
  unique

#creates directory for biolog raw data and sets it as working directory
dir.create("biolog_rawdata")
setwd('biolog_rawdata')

#downloads csv platereader files from google drive into biolog raw data folder
#this will open a window and ask for permissions to access your google drive
ls_tibble <- googledrive::drive_ls('https://drive.google.com/drive/u/1/folders/1N6NGWna7vIP8ht5pMW2QaedlS3ovzZ2n')
for (file_id in ls_tibble$id) {
  googledrive::drive_download(as_id(file_id))
}

#unzips zip files and deletes original zip files (just leaving csv files of plate reader data)
for (zipfile in list.files()){
  unzip(zipfile)
  unlink(zipfile)
}

#deletes a file that is not plate data
unlink('df_csv')

#creates an empty dataframe where averages will be appended in loop
absorbance_df_long=
  matrix(nrow=0, ncol=6) %>%
  as.data.frame() %>%
  set_colnames(c("Sample", "Assay", "Day", "Lid", "Substrate", "Absorbance"))

#creates a long format data frame of microplate data by appending absorbance values
#for each combination of day, sample, assay, and substrate.
for (plate_df in list.files()[grepl("ECO", list.files())]){
  
  #extracts day from file name
  day=strsplit(plate_df, "[[:punct:]]")[[1]][1]
  
  #extracts sample from file name
  sample=strsplit(plate_df, "[[:punct:]]")[[1]][3]
  
  #extracts assay (PM8 vs. ecoplate) from file name
  assay=strsplit(plate_df, "[[:punct:]]")[[1]][2]
  
  #creates variable for whether lid was on or off
  lid=case_when(strsplit(plate_df, "[[:punct:]]")[[1]][4]=="LIDON" ~ "LIDON",
                  strsplit(plate_df, "[[:punct:]]")[[1]][4]=="csv" ~ "LIDOFF")
  
  #reads in microplate df for each plate, dropping first column that includes plate letters
  #result is an 8 x 12 dataframe
  temp_df = read.csv(plate_df)[,-1] %>% as.matrix
  
  #loops through all substrates, calculates mean absorbance and appends to dataframe
  for (i in substrate_vector){
    #takes absorbance values for each replicate
    Abs1=temp_df[which(substrate_template[,1:4]==i)]
    Abs2=temp_df[which(substrate_template[,5:8]==i)]
    Abs3=temp_df[which(substrate_template[,9:12]==i)]
    
    #averages absorbance for 3 replicates
    avg_abs=mean(Abs1, Abs2, Abs3)
    
    newline <- data.frame(t(c(sample, assay, day, lid, i, avg_abs)))
    colnames(newline) <- colnames(absorbance_df_long)

    absorbance_df_long =
      rbind(absorbance_df_long, newline)
  }
}

#examine dataframe
head(absorbance_df_long)

#plots mean absorbance over for all substrates, 
#also creating a variable for which tree it came from
#the tree code "B" is for the blank
#note that the data has not been cleaned and there is still an effect of the
#lids for day 1 which
absorbance_df_long %>%
  mutate(Absorbance=as.numeric(Absorbance),
         Tree=str_extract(Sample, "^.")) %>%
  ggplot(aes(x=Day, y=Absorbance, group=Sample, color=Tree)) +
  geom_line() +
  facet_wrap(~Substrate, nrow=4)
