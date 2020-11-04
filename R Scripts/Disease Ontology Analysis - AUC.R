# Necessary libraries
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
library(tibble)
library(metaseqR)
library(conflicted)

# Set preferred packages
conflict_prefer("paste", "base")
conflict_prefer("desc", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("arrange", "dplyr")

# Load in datasets
load(here("Data", "Allen", "MTG_matrix_scaled.Rdata"))
load(here("Data", "export_data", "He_DS1_Human_averaged.Rdata"))
load(here("Data", "export_data", "Maynard_dataset_average.Rdata"))

###############################################################################

# Functions to use in Analysis #

# Adapted from Leon French - code to create gene sets from tmod using DisGeNet Disease Ontology annotations
create_geneBackground <- function(source_data) {
  geneBackground <- source_data %>%
    dplyr::select(gene_symbol) %>%
    pull ()
  names(geneBackground) <- geneBackground
  geneBackground <- sort(unique(names(as.list(geneBackground))))
}

geneSets_tmod <- function(dataset_geneBackground) {
  disgenet <- read_tsv(paste0(here("Data", "genelists", "DisGeNet", "curated_gene_disease_associations.tsv"))) 
  disgenet %<>% dplyr::select(symbol = geneSymbol, name = diseaseName, ID = diseaseId)
  disgenet %<>% filter(symbol %in% geneBackground)
  geneLists <- group_by(disgenet, ID) %>% 
    dplyr::summarise(name = paste(unique(name), collapse = ","), genes = unique(list(symbol)), size = dplyr::n()) %>% 
    filter(size >= minGOgroupSize & size <= maxGOgroupSize) 
  namedLists <- geneLists$genes
  names(namedLists) <- geneLists$ID
  idToName <- data.frame(ID = geneLists$ID, Title = geneLists$name)
  geneSets <- makeTmod(modules = idToName, modules2genes = namedLists)
  geneSetsTable <- tmod2DataFrame(geneSets)
}

# Function to generate an ordered list of genes arranged by descending order in 
# terms of gene expression per layer
generate_AUC_genelist <- function(data, layer) {
  layer_label <- layer
  data %<>%
    select(gene_symbol:Layer_6)
  sorted_data <- data[order(data[,layer_label], decreasing = TRUE),]
  genelist <- sorted_data %>% pull(gene_symbol)
  return(genelist)
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

# Widen scRNA-seq data by cell type
widen_scRNA_seq <- function(dataset) {
  wide_data <- dataset %>%
    mutate(cortical_layer_label = gsub("L", "Layer_", cortical_layer_label)) %>%
    ungroup() %>%
    pivot_wider(
      names_from = cortical_layer_label,
      values_from = mean_expression_scaled
    ) %>%
    group_by(class_label) %>%
    nest()
  return(wide_data)
}

# Function to generate a dataframe of sorted gene lists by layer for AUC
AUC_genelist_fn <- function(dataset) {
  # List of layers to create AUC genelists and tables for
  layer_list <- c("Layer_1", "Layer_2", "Layer_3", "Layer_4", "Layer_5", "Layer_6")
  
  AUC_genelist <- list()
  for (i in layer_list) {
    AUC_genelist[[i]] <- generate_AUC_genelist(dataset, i)
  }
  AUC_genelist_df <- as.data.frame(do.call(rbind, AUC_genelist)) %>%
    t() %>% as.data.frame()
}

# Function to generate list of tables of AUC values using sorted gene lists
AUC_table_list_fn <- function(AUC_genelist_df) {
  AUC_table_list <- list()
  for (i in 1:6) {
    AUC_table_list[[i]] <- AUC_table(AUC_genelist_df[,i])
  }
  return(AUC_table_list)
}

#Function to combine p-values across datasets joined by the same disease ontology term
combine_pval_AUC_table <- function(Maynard_table, He_table, MTG_table = NULL) {
  
  if (is.null(MTG_table)) {
    combined_table <- inner_join(Maynard_table, He_table, by = "MainTitle",
                                 suffix = c("_Maynard", "_He")) %>%
      column_to_rownames(var = "MainTitle") %>%
      dplyr::select(contains("P.value"), contains("AUC")) %>%
      filter(P.Value_He <= 1) %>%
      filter(P.Value_Maynard <= 1) 
  } else {
    combined_table <- inner_join(Maynard_table, He_table, by = "MainTitle",
                                 suffix = c("_Maynard", "_He"))
    combined_table <- inner_join(combined_table, MTG_table, by = "MainTitle") %>%
      column_to_rownames(var = "MainTitle") %>%
      rename(P.Value_Allen = P.Value, AUC_Allen = AUC, N1_Allen = N1, adj.P.Val_Allen = adj.P.Val, ID_Allen = ID, aspect_Allen = aspect,
             otherNames_Allen = otherNames, rank_Allen = rank) %>%
      dplyr::select(contains("P.value"), contains("AUC")) %>%
      filter(P.Value_He <= 1) %>%
      filter(P.Value_Maynard <= 1) 
  }
  
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

combined_AUC_table_fn <- function(Maynard_AUC_table, He_AUC_table, Allen_AUC_table, layer_range = c(1:6)) {
  combined_AUC_table_list <- list()
  for (i in layer_range) {
    combined_AUC_table_list[[i]] <- combine_pval_AUC_table(Maynard_AUC_table[[i]], He_AUC_table[[i]], Allen_AUC_table[[i]])
  }
  return(combined_AUC_table_list)
}

################################################################################

# Analysis script #

# Input parameters:
maxGOgroupSize <- 200
minGOgroupSize <- 10
goSource <- 'org.Hs.eg'

# Set up Maynard data
Maynard_geneBackground <- create_geneBackground(Maynard_dataset_average)
Maynard_geneSets <-  geneSets_tmod(Maynard_geneBackground)
# Create dataframe of genelist
Maynard_AUC_genelist_df <- AUC_genelist_fn(Maynard_dataset_average)
# Create tables of AUC values
Maynard_AUC_table_list <- AUC_table_list_fn(Maynard_AUC_genelist_df)

# Set up He data
He_geneBackground <- create_geneBackground(He_DS1_Human_averaged)
He_geneSets <-  geneSets_tmod(He_geneBackground)
# Create dataframe of genelist
He_AUC_genelist_df <- AUC_genelist_fn(He_DS1_Human_averaged)
# Create tables of AUC values
He_AUC_table_list <- AUC_table_list_fn(He_AUC_genelist_df)


# Widen and separate scRNA-seq data by cell type
MTG_matrix_wide <- widen_scRNA_seq(MTG_matrix_scaled)

# Create dataframes of genelists for each dataset
MTG_GABA_AUC_genelist_df <- AUC_genelist_fn(MTG_matrix_wide$data[[1]])
MTG_GLUT_AUC_genelist_df <- AUC_genelist_fn(MTG_matrix_wide$data[[2]])
MTG_NONN_AUC_genelist_df <- AUC_genelist_fn(MTG_matrix_wide$data[[3]])

# Create tables of AUC values
MTG_GABA_AUC_table_list <- AUC_table_list_fn(MTG_GABA_AUC_genelist_df)
MTG_GLUT_AUC_table_list <- AUC_table_list_fn(MTG_GLUT_AUC_genelist_df)
MTG_NONN_AUC_table_list <- AUC_table_list_fn(MTG_NONN_AUC_genelist_df)


#Combined tables of p-values for both He and Maynard data
combined_GABA_AUC_table_list <- combined_AUC_table_fn(Maynard_AUC_table_list, He_AUC_table_list, MTG_GABA_AUC_table_list)
combined_GLUT_AUC_table_list <- combined_AUC_table_fn(Maynard_AUC_table_list, He_AUC_table_list, MTG_GLUT_AUC_table_list, c(1:5)) # Error in fisher.method for layer 6
combined_NONN_AUC_table_list <- combined_AUC_table_fn(Maynard_AUC_table_list, He_AUC_table_list, MTG_NONN_AUC_table_list)


#Summarize top 20 disease ontology terms over an AUC value of 0.5 for both datasets, and surviving correction
print(head(filter(combined_GABA_AUC_table_list[[2]], (( (AUC_Maynard > 0.5) | (AUC_He) > 0.5 & (AUC_Allen > 0.5)) & adj_combined_p_value < 0.05)), n=20)) #insert combined tables


#Summarize bottom 20 disease ontology terms over an AUC value of 0.5 for both datasets, and surviving correction
print(head(filter(combined_NONN_AUC_table_list[[2]], ((AUC_Maynard < 0.5 | AUC_He < 0.5 | (AUC_Allen < 0.5) ) & adj_combined_p_value < 0.05/3)), n=20)) #insert combined tables


# Save .txt file of genes in "Mood Disorders" for gene ontology
mood_disorders_genelist <- geneSetsTable %>% filter(Title == "Mood Disorders") %>% pull(feature_id)
mood_disorder_genelist <- unlist(strsplit(mood_disorder_genelist, ","))
write.table(mood_disorder_genelist, file = here("Data", "genelists", "DO_mood.disorders.txt"), row.names = F, col.names = F, quote = F)

# Save .txt file of genes in "Epileptic encephalopathy" for gene ontology
epileptic_encephalopathy_genelist <- geneSetsTable %>% filter(Title == "Epileptic encephalopathy") %>% pull(feature_id)
epileptic_encephalopathy_genelist <- unlist(strsplit(epileptic_encephalopathy_genelist, ","))
write.table(epileptic_encephalopathy_genelist, file = here("Data", "genelists", "DO_epileptic.encephalopathy.txt"), row.names = F, col.names = F, quote = F)

# Save .txt file of genes in "Cocaine-related Disorders" for gene ontology
cocaine_related_disorders_genelist <- geneSetsTable %>% filter(Title == "Cocaine-Related Disorders") %>% pull(feature_id)
cocaine_related_disorders_genelist <- unlist(strsplit(cocaine_related_disorders_genelist, ","))
write.table(cocaine_related_disorders_genelist, file = here("Data", "genelists", "DO_cocaine.related.disorders.txt"), row.names = F, col.names = F, quote = F)

