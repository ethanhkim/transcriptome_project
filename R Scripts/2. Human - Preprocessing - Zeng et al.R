#Adapted from https://github.com/jritch/mri_transcriptomics/blob/master/R%20Code/RunSingleGO.AUROC.Analysis.R

#Load required libraries

library(tidyverse)
library(readxl)
library(dplyr)
library(magrittr)
library(HGNChelper)
library(here)
library(mygene)

Zeng_Path <- here("Data", "Zeng et al", "Table S2.xlsx")

Zeng_dataset <- read_xlsx(path = Zeng_Path, sheet = "Final1000New", skip=1)

Zeng_dataset %<>% 
  dplyr::select("Gene symbol", "Entrez Gene ID", "Cortical marker (human)", "Level...20", "Pattern...21",
                "Pattern...23", "Pattern...25") %>% 
  dplyr::rename(gene_symbol = "Gene symbol", entrez_id = "Entrez Gene ID", marker_annotation = "Cortical marker (human)", 
                expression_level = "Level...20", V1_pattern = "Pattern...21", V2_pattern = "Pattern...23", Temporal_pattern = "Pattern...25")

#Separate by layer marker
Zeng_dataset_expanded <- Zeng_dataset %>% 
  mutate(marker_annotation = strsplit(as.character(marker_annotation), "[/+]|( or )")) %>% 
  unnest(marker_annotation) %>%
  as_tibble()

#Assimilate layers
Zeng_dataset_expanded %<>% 
  mutate(marker_annotation = gsub("layer( )?","", marker_annotation)) %<>% 
  mutate(marker_annotation = gsub("[?]","", marker_annotation)) %<>% 
  mutate(marker_annotation = gsub("4c","4", marker_annotation)) %<>% 
  mutate(marker_annotation = gsub("5a","5", marker_annotation)) %<>% 
  mutate(marker_annotation = gsub("6b","6", marker_annotation)) %<>% 
  mutate(marker_annotation = gsub("([0-6])","layer \\1", marker_annotation)) %<>% 
  mutate(marker_annotation = gsub("VEC","vascular endothelial cell", marker_annotation)) %>%
  mutate(combined_annotation = paste(V1_pattern, V2_pattern, Temporal_pattern)) %>%
  mutate(combined_annotation = gsub("ubiquitous", "", combined_annotation)) %>%
  mutate(combined_annotation = gsub("layers", "", combined_annotation)) %>%
  mutate(combined_annotation = gsub("NA", "", combined_annotation)) %>%
  mutate(combined_annotation = gsub("sparse to none", "", combined_annotation)) %>%
  mutate(combined_annotation = gsub("ubiquitous", "", combined_annotation)) %>%
  mutate(combined_annotation = gsub("scattered", "", combined_annotation)) %>%
  mutate(combined_annotation = gsub("\\+", "", combined_annotation)) %>%
  mutate(combined_annotation = gsub("sparse", "", combined_annotation)) %>%
  mutate(combined_annotation = gsub("\\(", "", combined_annotation)) %>%
  mutate(combined_annotation = gsub("\\)", "", combined_annotation)) %>%
  mutate(combined_annotation = gsub("\\?", "", combined_annotation)) %>%
  mutate(combined_annotation = gsub("layer", "", combined_annotation)) %>%
  mutate(combined_annotation = gsub("\\,", " ", combined_annotation)) %>%
  mutate(combined_annotation = gsub("\\/", "", combined_annotation)) %>%
  mutate(combined_annotation = gsub("ubiquitous", "", combined_annotation)) %>%

Zeng_dataset_cleaned <- Zeng_dataset_expanded %>% 
  dplyr::filter(marker_annotation  != 'others' & marker_annotation != 'laminar' | is.na(NA))


#Check if gene_symbol list is updated
updated_symbols <- mygene::getGenes(Zeng_dataset_cleaned$entrez_id) %>% as.data.frame() %>% dplyr::select(entrez_id = entrezgene, updated_symbol = symbol)


#Update gene_symbol list using mygene
Zeng_dataset_updated <- Zeng_dataset_cleaned
Zeng_dataset_updated %<>% mutate(entrez_id = as.character(entrez_id))
#NCBI ID change for one gene 10571 -> 11039
Zeng_dataset_updated %<>% mutate(entrez_id = if_else(entrez_id == "10571", "11039", entrez_id))
Zeng_dataset_updated <- left_join(Zeng_dataset_updated , updated_symbols)
Zeng_dataset_updated %<>% mutate(gene_symbol = if_else(!is.na(updated_symbol), updated_symbol, gene_symbol))
Zeng_dataset_updated %<>% dplyr::select(-updated_symbol)


#Make into long format
Zeng_dataset_updated$entrez_id <- as.character(Zeng_dataset_updated$entrez_id)

Zeng_dataset_long <- Zeng_dataset_updated %>%
  pivot_longer(
    cols = V1_pattern:Temporal_pattern,
    names_to = "region",
    values_to = "original_annotation"
  ) %>%
  mutate(region = gsub("V1_pattern", "V1", region)) %>%
  mutate(region = gsub("V2_pattern", "V2", region)) %>%
  mutate(region = gsub("Temporal_pattern", "Temporal", region)) %>%
  dplyr::select("gene_symbol", "entrez_id", "region", "expression_level", "marker_annotation", "original_annotation")

#Summarize results
Zeng_dataset_long %>% 
  group_by(marker_annotation) %>% 
  summarise(n())

write_csv(Zeng_dataset_updated, './R Scripts/export_data/Cleaned_Zeng_dataset.csv')
save(Zeng_dataset_updated, file = './R Scripts/export_data/Zeng_dataset_updated.Rdata')

write_csv(Zeng_dataset_long, './R Scripts/export_data/Zeng_dataset_long.csv')
save(Zeng_dataset_long, file = './R Scripts/export_data/Zeng_dataset_long.Rdata')
