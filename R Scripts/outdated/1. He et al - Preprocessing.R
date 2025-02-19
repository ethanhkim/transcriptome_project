#Preprocessing He et al data to change Ensembl ID to Gene ID

#Load required libraries

library(tidyverse)
library(org.Hs.eg.db)
library(here)
library(readxl)
library(dplyr)
library(HGNChelper)
library(conflicted)

# Set conflicts
conflict_prefer('filter', 'dplyr')
conflict_prefer('rename', 'dplyr')
conflict_prefer('select', 'dplyr')

## DS1 ##
He_DS1_Human <- read.table(here("Data", "raw_data", "He et al", 
                                "tab.average_RPKM_section_human.tsv"), 
                           header = T)

#Truncate the version numbers on ENSEMBL ID
He_DS1_Human$ENSEMBL_ID <- gsub("\\.\\d+$", "", He_DS1_Human$ENSEMBL_ID)

#Change ENSEMBL ID to Gene symbol
He_DS1_Human$gene_symbol <- mapIds(org.Hs.eg.db, 
                                   keys = He_DS1_Human$ENSEMBL_ID, 
                                   keytype = "ENSEMBL", column="SYMBOL")

#Check if genes are updated
check_HeList <- checkGeneSymbols(He_DS1_Human$gene_symbol, 
                                 unmapped.as.na = TRUE, map = NULL, 
                                 species = "human")

#Update gene_symbol list
He_DS1_Human <- He_DS1_Human %>%
  add_column(updated_gene_symbol = check_HeList$Suggested.Symbol) %>%
  dplyr::select(-"gene_symbol") %>%
  dplyr::rename(gene_symbol = "updated_gene_symbol")

#Rename
He_DS1_Human <- He_DS1_Human %>%
  rename_at(vars(contains('X')), funs(sub('X', 'S', .))) %>%
  rename(Ensembl_ID = 'ENSEMBL_ID')

#Remove unnecessary columns, rows and reorder so gene_symbol is at the front
He_DS1_Human <- He_DS1_Human %>%
  select(-"Ensembl_ID") %>%
  select("gene_symbol", everything()) %>%
  filter(!is.na(gene_symbol)) %>%
  #Remove duplicates
  distinct(gene_symbol, .keep_all = TRUE)

# Clean up environment
rm(check_HeList)

## DS2 not used in analysis ##

#DS2_Human <- read.table(here("Data", "He et al", "DS2_sectionRPKM_human.tsv"), stringsAsFactors = FALSE)
#DS2_Human$V1 <- gsub("\\.\\d+$", "", DS2_Human$V1)
#DS2_Human$gene_symbol <- mapIds(org.Hs.eg.db, keys = DS2_Human$V1, keytype = "ENSEMBL", column="SYMBOL")

#DS2_Human <- DS2_Human %>%
  #rename(Ensembl_ID = "V1",
         #V1 = "V2", V2 = "V3", V3 = "V4", V4 = "V5", V5 = "V6", V6 = "V7", V7 = "V8",
         #V8 = "V9", V9 = "V10", V10 = "V11") 

#DS2_Human <- DS2_Human %>%
  #dplyr::select(-"Ensembl_ID") %>%
  #dplyr::select("gene_symbol", everything()) %>%
  #filter(!is.na(gene_symbol)) %>%
  #distinct(gene_symbol, .keep_all = TRUE)
