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
library(rvest)
library(httr)

# Set preferred packages
conflict_prefer("paste", "base")
conflict_prefer("desc", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("arrange", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("strsplit", "base")
conflict_prefer("content", "httr")

# Load in datasets
load(here("Data", "processed_data", "Allen_logCPM_dataset.Rdata"))
load(here("Data", "processed_data", "Allen_logCPM_filtered_dataset.Rdata"))
load(here("Data", "processed_data", "Allen_downsampled_logCPM_dataset.Rdata"))
load(here("Data", "processed_data", "Allen_downsampled_logCPM_filtered_dataset.Rdata"))
load(here("Data", "processed_data", "Maynard_logCPM_dataset.Rdata"))
load(here("Data", "processed_data", "Maynard_logCPM_filtered_dataset.Rdata"))
load(here("Data", "processed_data", "He_DS1_logCPM_dataset.Rdata"))
load(here("Data", "processed_data", "He_DS1_logCPM_filtered_dataset.Rdata"))

# Source functions
source(here('R Scripts','Layer-specific analysis', 'Gene and Disease Ontology Analysis functions.R'))

# Analysis script ----

# Input parameters:
maxGOgroupSize <- 200
minGOgroupSize <- 10
goSource <- 'org.Hs.eg'
ontology_type <- "gene" # input either "gene" or "disease" depending on if wanting to analyze gene/disease ontology

# Create list of tables of AUC values for DO groups for Maynard
Maynard_AUC_table_list <- dataset_AUC_summary(Maynard_logCPM_filtered_dataset, ontology_type)
# Create list of tables of AUC values for DO groups for He
He_AUC_table_list <- dataset_AUC_summary(He_DS1_logCPM_filtered_dataset, ontology_type)

# Widen and separate scRNA-seq data by cell type
Allen_matrix_wide <- widen_scRNA_seq(Allen_downsampled_logCPM_filtered_dataset)

# Create dataframes of genelists for each single cell dataset
Allen_GABA_AUC_table_list <- dataset_AUC_summary(Allen_matrix_wide$data[[1]], ontology_type)
Allen_GLUT_AUC_table_list <- dataset_AUC_summary(Allen_matrix_wide$data[[2]], ontology_type)
Allen_NONN_AUC_table_list <- dataset_AUC_summary(Allen_matrix_wide$data[[3]], ontology_type)

# Combined tables of p-values for both He, Maynard and Allen data
combined_GABA_AUC_table_list <- combined_AUC_table_fn(Maynard_AUC_table_list, He_AUC_table_list, Allen_GABA_AUC_table_list)
combined_GLUT_AUC_table_list <- combined_AUC_table_fn(Maynard_AUC_table_list, He_AUC_table_list, Allen_GLUT_AUC_table_list)
combined_NONN_AUC_table_list <- combined_AUC_table_fn(Maynard_AUC_table_list, He_AUC_table_list, Allen_NONN_AUC_table_list)

# Ontology tables
GABA_ontology <- extract_term_and_pvalue(combined_GABA_AUC_table_list, 20) %>%
  add_column(cell_type = "GABA")
GLUT_ontology <- extract_term_and_pvalue(combined_GLUT_AUC_table_list, 20) %>%
  add_column(cell_type = "GLUT")
NONN_ontology <- extract_term_and_pvalue(combined_NONN_AUC_table_list, 20) %>%
  add_column(cell_type = "NONN")
# Create common layer ontology term table
common_ontology_layer <- rbind(GABA_ontology, GLUT_ontology, NONN_ontology) %>%
  filter(rank <= 20) %>%
  group_by(ID) %>%
  filter(n()>1) %>%
  arrange(Ontology_Term) %>%
  select(rank, Ontology_Term, ID, mean_AUC, AUC_Maynard, AUC_He, 
         AUC_Allen, combined_p_value, direction, layer, cell_type)

# For each layer, print the enriched terms within the top 20 rank
for (layer_number in c(1, 2, 3, 4, 5, 6)) {
  print(common_ontology_layer %>%
          filter(layer == layer_number & direction == "positive" & rank <= 20) %>%
          group_by(ID) %>%
          filter(n()>1) %>%
          arrange(Ontology_Term))
}
# For each layer, print the depleted terms within the top 20 rank
for (layer_number in c(1, 2, 3, 4, 5, 6)) {
  print(common_ontology_layer %>%
          filter(layer == layer_number & direction == "negative" & rank <= 20) %>%
          group_by(ID) %>%
          filter(n()>1) %>%
          arrange(Ontology_Term))
}

View(common_ontology_layer %>%
       filter(layer == 5 & direction == "negative" & rank <= 20) %>%
       group_by(ID) %>%
       filter(n()>1) %>%
       arrange(Ontology_Term))

write.csv(common_ontology_layer, file = here("Data", "ontology_tables", "CPM_nonDownsampled_gene_ontology.csv"))



# ReviGo
#for (layer_number in c(1, 2, 3, 4, 5, 6)) {
#  run_revigo(common_ontology_layer, layer_number, "positive")
#}
#run_revigo(GABA_ontology, layer_number = 1, "positive")

#run_revigo()




#Summarize top and bottom 20 disease ontology terms over an AUC value of 0.5 for both datasets, and surviving correction ----
#print_ontology_terms(combined_GABA_AUC_table_list, 3, "positive")
#print_ontology_terms(combined_GLUT_AUC_table_list, 3, "positive")
#print_ontology_terms(combined_NONN_AUC_table_list, 3, "positive")

#Summarize top and bottom 20 gene ontology terms over an AUC value of 0.5 for both datasets, and surviving correction ----
#print_ontology_terms(combined_GABA_AUC_table_list, 3, "negative")
#print_ontology_terms(combined_GLUT_AUC_table_list, 3, "negative")
#print_ontology_terms(combined_NONN_AUC_table_list, 3, "negative")

#Summarize top and bottom 20 disease ontology terms over an AUC value of 0.5 for both datasets, and surviving correction
#print(head(filter(GO_GABA_AUC_table[[3]], (( (AUC_Maynard > 0.5) & (AUC_He > 0.5) & (AUC_Allen > 0.5)) & adj_combined_p_value < 0.05 & P.Value_Allen < 0.05)), n=20)) #insert combined tables
#Summarize bottom 20 disease ontology terms over an AUC value of 0.5 for both datasets, and surviving correction
#print(head(filter(GO_NONN_AUC_table[[2]], (( (AUC_Maynard < 0.5) & (AUC_He < 0.5) & (AUC_Allen < 0.5)) & adj_combined_p_value < 0.05 & P.Value_Allen < 0.05)), n=20))

# Save .txt file of genes in "Mood Disorders" for gene ontology
#mood_disorders_genelist <- geneSetsTable %>% filter(Title == "Mood Disorders") %>% pull(feature_id)
#mood_disorder_genelist <- unlist(strsplit(mood_disorder_genelist, ","))
#write.table(mood_disorder_genelist, file = here("Data", "genelists", "DO_mood.disorders.txt"), row.names = F, col.names = F, quote = F)

# Save .txt file of genes in "Epileptic encephalopathy" for gene ontology
#epileptic_encephalopathy_genelist <- geneSetsTable %>% filter(Title == "Epileptic encephalopathy") %>% pull(feature_id)
#epileptic_encephalopathy_genelist <- unlist(strsplit(epileptic_encephalopathy_genelist, ","))
#write.table(epileptic_encephalopathy_genelist, file = here("Data", "genelists", "DO_epileptic.encephalopathy.txt"), row.names = F, col.names = F, quote = F)

# Save .txt file of genes in "Cocaine-related Disorders" for gene ontology
#cocaine_related_disorders_genelist <- geneSetsTable %>% filter(Title == "Cocaine-Related Disorders") %>% pull(feature_id)
#cocaine_related_disorders_genelist <- unlist(strsplit(cocaine_related_disorders_genelist, ","))
#write.table(cocaine_related_disorders_genelist, file = here("Data", "genelists", "DO_cocaine.related.disorders.txt"), row.names = F, col.names = F, quote = F)


## ReviGO Analyses ----
ReviGO_GABA_list_up <- list()
for (i in 1:6) {
  ReviGO_GABA_list_up[[i]] <- run_revigo(ReviGO_GABA, i, "positive")
}
ReviGO_GABA_list_down <- list()
for (i in 1:6) {
  ReviGO_GABA_list_down[[i]] <- run_revigo(ReviGO_GABA, i, "negative")
}

ReviGO_GABA_up <- do.call(rbind, ReviGO_GABA_list_up)
ReviGO_GABA_down <- do.call(rbind, ReviGO_GABA_list_down)
ReviGO_GABA_complete <- rbind(ReviGO_GABA_up, ReviGO_GABA_down)

ReviGO_GLUT_list_up <- list()
for (i in 1:6) {
  ReviGO_GLUT_list_up[[i]] <- run_revigo(ReviGO_GLUT, i, "positive")
}
ReviGO_GLUT_list_down <- list()
for (i in 1:6) {
  ReviGO_GLUT_list_down[[i]] <- run_revigo(ReviGO_GLUT, i, "negative")
}

ReviGO_GLUT_up <- do.call(rbind, ReviGO_GLUT_list_up)
ReviGO_GLUT_down <- do.call(rbind, ReviGO_GLUT_list_down)
ReviGO_GLUT_complete <- rbind(ReviGO_GLUT_up, ReviGO_GLUT_down)

ReviGO_NONN_list_up <- list()
for (i in 1:6) {
  ReviGO_NONN_list_up[[i]] <- run_revigo(ReviGO_NONN, i, "positive")
}
ReviGO_NONN_list_down <- list()
for (i in 1:6) {
  ReviGO_NONN_list_down[[i]] <- run_revigo(ReviGO_NONN, i, "negative")
}

ReviGO_NONN_up <- do.call(rbind, ReviGO_NONN_list_up)
ReviGO_NONN_down <- do.call(rbind, ReviGO_NONN_list_down)
ReviGO_NONN_complete <- rbind(ReviGO_NONN_up, ReviGO_NONN_down)

ReviGO_GABA[1]


# SRP9 gene - look at allen brain mouse atlas and Zeng images 