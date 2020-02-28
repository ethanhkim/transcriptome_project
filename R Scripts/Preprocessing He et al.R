#Preprocessing He et al data to change Ensembl ID to Gene ID

#Load required libraries

library(tidyverse)
library(org.Hs.eg.db)
library(org.Mmu.eg.db)

DS2_Human <- read.table("./Data/He et al/DS2_sectionRPKM_human.tsv")
DS2_Macaque <- read.table("./Data/He et al/Ds2_sectionRPKM_macaque.tsv")
DS1_Human <- read.table("./Data/He et al/tab.average_RPKM_section_human.tsv")
DS1_Macaque <- read.table("./Data/He et al/tab.average_RPKM_section_macaque.tsv")
DS1_Chimp <- read.table("./Data/He et al/tab.average_RPKM_section_chimp.tsv")

#Provided Ensembl ID are all the same; refer to error here:

#Truncate the version numbers on ENSEMBL ID
DS2_Human$V1 <- gsub("\\.\\d+$", "", DS2_Human$V1)
DS2_Macaque$V1 <- gsub("\\.\\d+$", "", DS2_Macaque$V1)
DS1_Human$V1 <- gsub("\\.\\d+$", "", DS1_Human$V1)
DS1_Macaque$V1 <- gsub("\\.\\d+$", "", DS1_Macaque$V1)
DS1_Chimp$V1 <- gsub("\\.\\d+$", "", DS1_Chimp$V1)

#Change ENSEMBL ID to Gene ID
#Note: needs for loop
DS2_Human$gene_ID <- mapIds(org.Hs.eg.db, keys = DS2_Human$V1, keytype = "ENSEMBL", column="SYMBOL")
DS2_Macaque$geneID <- mapIds(org.Hs.eg.db, keys = DS2_Macaque$V1, keytype = "ENSEMBL", column="SYMBOL")
DS1_Human$gene_ID <- mapIds(org.Hs.eg.db, keys = DS1_Human$V1, keytype = "ENSEMBL", column="SYMBOL")
DS1_Macaque$gene_ID <- mapIds(org.Hs.eg.db, keys = DS1_Macaque$V1, keytype = "ENSEMBL", column="SYMBOL")
DS1_Chimp$gene_ID <- mapIds(org.Hs.eg.db, keys = DS1_Chimp$V1, keytype = "ENSEMBL", column="SYMBOL")

#Reorder
DS1_Human <- DS1_Human[, c(1, 19, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,13,14,15,16,17,18)]
DS1_Chimp <- DS1_Chimp[, c(1, 19, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,13,14,15,16,17,18)]
DS1_Macaque <- DS1_Macaque[, c(1, 19, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,13,14,15,16,17,18)]
DS2_Human <- DS2_Human[, c(1, 12, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)]
DS2_Macaque <- DS2_Macaque[, c(1, 12, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)]
