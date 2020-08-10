library(cowplot)
library(readr)
library(ggplot2)
library(tidyr)
library(magrittr)
library(AnnotationDbi)
library(annotate)
library(GO.db)
library(tmod)
library(org.Hs.eg.db) 
library(here)
library(dplyr)
library(metaseqR)
library(tidyverse)
library(data.table)

#Input parameters:
autismGeneListsFolder <- here("Data", "genelists", "autism genelists")

# Adapted from Leon French - code to create gene sets from tmod using curated autism gene lists
loadFileSets <- function(prefix) {
  tmodNames <- data.frame()
  modules2genes <- list()
  
  for(geneListFilename in list.files(autismGeneListsFolder, pattern = paste0(prefix,".*txt"), full.names = T)) {
    print(geneListFilename)
    genesOfInterest <- read.csv(geneListFilename,header=F,stringsAsFactors = F)
    shortName <- gsub(".txt","",gsub(paste0(".*/"),"", geneListFilename))
    shortName <- gsub(paste0(prefix,"."),"",shortName)
    
    genesOfInterest$term <- shortName
    
    modules2genes[shortName] <- list(genesOfInterest$V1)
    
    tmodNames <- rbind(tmodNames, data.frame(ID=shortName, Title = shortName))
  }
  geneSets <- makeTmod(modules = tmodNames, modules2genes = modules2genes)
  geneSets
}


#Function to generate an ordered list of genes arranged by descending order in terms of gene expression per layer
generate_AUC_genelist <- function(data, layer) {
  data %>%
    dplyr::select(gene_symbol:Layer_6) %>%
    dplyr::arrange(desc({{ layer }})) %>%
    dplyr::select(gene_symbol) %>%
    pull()
}

#Function to create table of AUC and p-values associated with disease ontology terms
AUC_table <- function(AUC_genelist) {
  #Creates a top and bottom 20 list of ontology definitions using AUC
  AUC_table <- tbl_df(tmodUtest(c(AUC_genelist), mset=geneSets, qval = 1.01, filter = T))
  AUC_table %<>% rowwise() %>% mutate(P.Value = P.Value * 2) %>% ungroup() %>% mutate(adj.P.Val=p.adjust(P.Value, method="fdr")) #tmod runs one-sided tests
  AUC_table %<>% rowwise() %>% mutate(aspect = Ontology(ID)) #add the source ontology (could be filtered for just biological process)
  
  #collapse genesets that have the exact same set of genes
  AUC_table %<>% rowwise() %>% mutate(genes = paste(sort(unlist(geneSets$MODULES2GENES[ID])), collapse = " "))
  AUC_table %<>% ungroup() %>% group_by(genes, N1) %>% arrange(Title) %>% 
    summarize(MainTitle = dplyr::first(Title),  ID=paste(ID, collapse=","), AUC = dplyr::first(AUC), 
              P.Value= dplyr::first(P.Value), aspect= dplyr::first(aspect), 
              otherNames = if_else(length(unique(Title)) > 1, paste(Title[2:length(Title)], collapse=", "), ""))
  AUC_table %<>% ungroup() %>% mutate(adj.P.Val=p.adjust(P.Value, method="fdr"))
  AUC_table %<>% dplyr::select(-genes)
  AUC_table %<>% arrange(P.Value)
  AUC_table$rank <- 1:nrow(AUC_table)
  AUC_table %<>% dplyr::select(MainTitle, N1, AUC, P.Value, adj.P.Val, everything()) %>% arrange(adj.P.Val)
}

#Function to combine p-values across datasets joined by the same disease ontology term
combine_pval_AUC_table <- function(Maynard_table, He_table) {
  
  combined_table <- inner_join(Maynard_table, He_table, by = "MainTitle",
                               suffix = c("_Maynard", "_He")) %>%
    column_to_rownames(var = "MainTitle") %>%
    dplyr::select(contains("P.value"), contains("AUC")) %>%
    filter(P.Value_He <= 1) %>%
    filter(P.Value_Maynard <= 1)
  
  p_value_fisher <- fisher.method(combined_table, method = c("fisher"))
  #adjust the p-value using fdr
  p_value_fisher$adj_p_fdr <- p.adjust(p_value_fisher$p.value, method = "fdr")
  
  combined_table$combined_p_value <- p_value_fisher$p.value
  combined_table$adj_combined_p_value <- p_value_fisher$adj_p_fdr
  
  combined_table <- combined_table %>%
    rownames_to_column(var = "Disease_Ontology_Term") %>%
    arrange(adj_combined_p_value)
  combined_table$rank <- 1:nrow(combined_table)
  return(combined_table)
}

# Create autism gene sets
autism_sets <- loadFileSets("Autism")

# Create Maynard AUC gene lists by layer
Maynard_AUC_genelist_Layer1 <- generate_AUC_genelist(Maynard_dataset_average, Layer_1)
Maynard_AUC_genelist_Layer2 <- generate_AUC_genelist(Maynard_dataset_average, Layer_2)
Maynard_AUC_genelist_Layer3 <- generate_AUC_genelist(Maynard_dataset_average, Layer_3)
Maynard_AUC_genelist_Layer4 <- generate_AUC_genelist(Maynard_dataset_average, Layer_4)
Maynard_AUC_genelist_Layer5 <- generate_AUC_genelist(Maynard_dataset_average, Layer_5)
Maynard_AUC_genelist_Layer6 <- generate_AUC_genelist(Maynard_dataset_average, Layer_6)

# Create He et al AUC gene lists by layer
He_AUC_genelist_Layer1 <- generate_AUC_genelist(He_DS1_Human_averaged, Layer_1)
He_AUC_genelist_Layer2 <- generate_AUC_genelist(He_DS1_Human_averaged, Layer_2)
He_AUC_genelist_Layer3 <- generate_AUC_genelist(He_DS1_Human_averaged, Layer_3)
He_AUC_genelist_Layer4 <- generate_AUC_genelist(He_DS1_Human_averaged, Layer_4)
He_AUC_genelist_Layer5 <- generate_AUC_genelist(He_DS1_Human_averaged, Layer_5)
He_AUC_genelist_Layer6 <- generate_AUC_genelist(He_DS1_Human_averaged, Layer_6)

#Maynard et al AUC tables
Maynard_AUC_table_Layer1 <- AUC_table(Maynard_AUC_genelist_Layer1)
Maynard_AUC_table_Layer2 <- AUC_table(Maynard_AUC_genelist_Layer2)
Maynard_AUC_table_Layer3 <- AUC_table(Maynard_AUC_genelist_Layer3)
Maynard_AUC_table_Layer4 <- AUC_table(Maynard_AUC_genelist_Layer4)
Maynard_AUC_table_Layer5 <- AUC_table(Maynard_AUC_genelist_Layer5)
Maynard_AUC_table_Layer6 <- AUC_table(Maynard_AUC_genelist_Layer6)

#He et al AUC tables
He_AUC_table_Layer1 <- AUC_table(He_AUC_genelist_Layer1)
He_AUC_table_Layer2 <- AUC_table(He_AUC_genelist_Layer2)
He_AUC_table_Layer3 <- AUC_table(He_AUC_genelist_Layer3)
He_AUC_table_Layer4 <- AUC_table(He_AUC_genelist_Layer4)
He_AUC_table_Layer5 <- AUC_table(He_AUC_genelist_Layer5)
He_AUC_table_Layer6 <- AUC_table(He_AUC_genelist_Layer6)

#Combined tables of p-values for both He and Maynard data
combined_table_Layer1 <- combine_pval_AUC_table(Maynard_AUC_table_Layer1, He_AUC_table_Layer1) 
combined_table_Layer2 <- combine_pval_AUC_table(Maynard_AUC_table_Layer2, He_AUC_table_Layer2) 
combined_table_Layer3 <- combine_pval_AUC_table(Maynard_AUC_table_Layer3, He_AUC_table_Layer3)
combined_table_Layer4 <- combine_pval_AUC_table(Maynard_AUC_table_Layer4, He_AUC_table_Layer4)
combined_table_Layer5 <- combine_pval_AUC_table(Maynard_AUC_table_Layer5, He_AUC_table_Layer5)
combined_table_Layer6 <- combine_pval_AUC_table(Maynard_AUC_table_Layer6, He_AUC_table_Layer6)