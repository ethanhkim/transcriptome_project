#Adapted from https://github.com/jritch/mri_transcriptomics/blob/master/R%20Code/RunSingleGO.AUROC.Analysis.R

#Load required libraries

library(tidyverse)
library(readxl)
library(dplyr)
library(magrittr)

zengPath <- "./Data/Zeng et al/Table S2.xlsx"

zengTable <- read_xlsx(path = zengPath, sheet = "Final1000New", skip=1)

zengTable %<>% 
  dplyr::select("Gene symbol", "Cortical marker (human)", "Level...20") %>% 
  dplyr::rename(Gene.symbol = "Gene symbol", Cortical.marker..human. = "Cortical marker (human)", Expression.level = "Level...20")

#Separate by layer marker
expandedZengTable <- zengTable %>% 
  mutate(Cortical.marker..human. = strsplit(as.character(Cortical.marker..human.), "[/+]|( or )")) %>% 
  unnest(Cortical.marker..human.)

as.data.frame(expandedZengTable)

#Assimilate layers
expandedZengTable %<>% 
  mutate(Cortical.marker..human. = gsub("layer( )?","", Cortical.marker..human.)) %<>% 
  mutate(Cortical.marker..human. = gsub("[?]","", Cortical.marker..human.)) %<>% 
  mutate(Cortical.marker..human. = gsub("4c","4", Cortical.marker..human.)) %<>% 
  mutate(Cortical.marker..human. = gsub("5a","5", Cortical.marker..human.)) %<>% 
  mutate(Cortical.marker..human. = gsub("6b","6", Cortical.marker..human.)) %<>% 
  mutate(Cortical.marker..human. = gsub("([0-6])","layer \\1", Cortical.marker..human.)) %<>% 
  mutate(Cortical.marker..human. = gsub("VEC","vascular endothelial cell", Cortical.marker..human.))

cleanedZengTable <- expandedZengTable %>% 
  dplyr::filter(Cortical.marker..human. != 'others' & Cortical.marker..human. != "laminar" | is.na(NA))

cleanedZengTable %>% 
  group_by(Cortical.marker..human.) %>% 
  summarise(n())
