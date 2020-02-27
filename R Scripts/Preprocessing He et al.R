#Preprocessing He et al

#Load required libraries

library(tidyverse)
library(org.Hs.eg.db)
library(biomaRt)
library(org.Mmu.eg.db)

DS2_Human <- read.table("./Data/He et al/DS2_sectionRPKM_human.tsv")
DS2_Macaque <- read.table("./Data/He et al/Ds2_sectionRPKM_macaque.tsv")
DS1_Human <- read.table("./Data/He et al/tab.average_RPKM_section_human.tsv")
DS1_Macaque <- read.table("./Data/He et al/tab.average_RPKM_section_macaque.tsv")
DS1_Chimp <- read.table("./Data/He et al/tab.average_RPKM_section_chimp.tsv")

#Extract Ensembl ID from row
DS1_Human_ensembl_id <- rownames(DS1_Human)
DS1_Human$ENSEMBL_ID <- DS1_Human_ensembl_id

DS1_Macaque_ensembl_id <- rownames(DS1_Macaque)
DS1_Macaque$ENSEMBL_ID <- DS1_Macaque_ensembl_id

#Truncate the version numbers on ENSEMBL ID
DS2_Human_RPKM_cleaned <- gsub("\\.\\d+$", "", DS2_Human$V1)
DS2_Macaque_cleaned <- gsub("\\.\\d+$", "", DS2_Macaque$V1)
DS1_Human$ENSEMBL_ID <- gsub("\\.\\d+$", "", DS1_Human$ENSEMBL_ID)
DS1_Macaque$ENSEMBL_ID <- gsub("\\.\\d+$", "", DS1_Macaque$ENSEMBL_ID)

#Merge cleaned column into main TSV
DS2_Human$V1 <- DS2_Human_RPKM_cleaned
DS2_Macaque$V1 <- DS2_Macaque_cleaned

#Change ENSEMBL ID to Gene ID for human
DS2_human_geneID <- mapIds(org.Hs.eg.db, keys = DS2_Human$V1, keytype = "ENSEMBL", column="SYMBOL")
DS2_macaque_geneID <- mapIds(org.Hs.eg.db, keys = DS2_Macaque$V1, keytype = "ENSEMBL", column="SYMBOL")
DS1_Human$gene_ID <- mapIds(org.Hs.eg.db, keys = DS1_Human$ENSEMBL_ID, keytype = "ENSEMBL", column="SYMBOL")
DS1_Macaque$gene_ID <- mapIds(org.Hs.eg.db, keys = DS1_Macaque$ENSEMBL_ID, keytype = "ENSEMBL", column="SYMBOL")

#Add ENSEMBL ID to table
DS2_Human$Gene_ID <- DS2_human_geneID
DS2_Macaque$Gene_ID <- DS2_macaque_geneID