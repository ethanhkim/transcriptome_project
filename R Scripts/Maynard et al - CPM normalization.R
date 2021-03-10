## Maynard et al CPM normalization ##

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("LieberInstitute/spatialLIBD")
library(spatialLIBD)
library(data.table)
library(org.Hs.eg.db)
library(dplyr)
library(magrittr)
library(stringr)
library(conflicted)
library(biomaRt)
library(HGNChelper)
library(edgeR)
library(here)

# Set conflicts #
conflict_prefer('select', 'dplyr')
conflict_prefer('filter', 'dplyr')
conflict_prefer('cpm', 'edgeR')
conflict_prefer('rename', 'dplyr')

# Fetch data from spatialLIBD
sce_layer <- fetch_data(type = "sce_layer")
# Get raw UMI counts per gene
Maynard_dataset <- sce_layer@assays@data@listData$counts %>%
  as_tibble(rownames=NA)
# Get gene Ensembl ID
Maynard_ensembl_list <- sce_layer@rowRanges@ranges@NAMES
# Convert Ensembl to HGNC, filter out non-viable symbols
Maynard_dataset$gene_symbol <- mapIds(org.Hs.eg.db, 
                                      keys = Maynard_ensembl_list, keytype = "ENSEMBL", column="SYMBOL")
Maynard_dataset %<>%
  select(gene_symbol, everything()) %>% filter(!is.na(gene_symbol)) %>%
  distinct(gene_symbol, .keep_all = TRUE)

# Select only individuals with all cortical layers (n = 2)
Maynard_dataset_subset <- Maynard_dataset %>%
  select("gene_symbol", contains('151507'), contains('151508'), contains('151509'),
         contains('151510'), contains('151673'), contains('151674'), 
         contains('151675'), contains('151676')) %>%
  column_to_rownames(var = "gene_symbol")

# Function to create columns of mean values across individuals for Maynard
create_mean_col <- function(Maynard_normalized_subset, layer_label) {
  df <- Maynard_normalized_subset
  layer <- layer_label
  
  df %<>% select(contains(layer)) %>%
    mutate(mean = rowMeans(.)) %>%
    pull(mean)
  
  return(df)
}

# Create list of mean columns
Maynard_mean_col_list <- list()
for (label in c("Layer1", "Layer2", "Layer3", "Layer4", "Layer5", "Layer6", "WM")) {
  Maynard_mean_col_list[[label]] <- create_mean_col(Maynard_dataset_subset, label)
}
# Bind list and convert to dataframe
Maynard_dataset_averaged <- as.data.frame(do.call(rbind, Maynard_mean_col_list)) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var = 'gene_symbol') %>%
  rename(Layer_1 = Layer1, Layer_2 = Layer2, Layer_3 = Layer3, Layer_4 = Layer4, 
         Layer_5 = Layer5, Layer_6 = Layer6)

# Normalize Maynard UMI counts with CPM, log = T
Maynard_logCPM_dataset <- Maynard_dataset_averaged %>%
  # Remove gene_symbol column for cpm()
  select(-gene_symbol) %>% cpm(log = T) %>%
  as.data.frame() %>%
  add_column(gene_symbol = Maynard_dataset_averaged$gene_symbol) %>%
  select(gene_symbol, everything())

# Normalize Maynard UMI counts with CPM, log = T
Maynard_CPM_dataset <- Maynard_dataset_averaged %>%
  # Remove gene_symbol column for cpm()
  select(-gene_symbol) %>% cpm() %>%
  as.data.frame() %>%
  add_column(gene_symbol = Maynard_dataset_averaged$gene_symbol) %>%
  select(gene_symbol, everything())

# Clean up workspace
rm(Maynard_dataset, Maynard_dataset_averaged, Maynard_dataset_subset,
   Maynard_mean_col_list, label, Maynard_ensembl_list,
   create_mean_col, sce_layer)

# Write normalized data as .Rdata
save(Maynard_logCPM_dataset, file = here("Data", "processed_data", 
                                         "Maynard_logCPM_dataset.Rdata"))
save(Maynard_CPM_dataset, file = here("Data", "processed_data", 
                                      "Maynard_CPM_dataset.Rdata"))

# Write normalized data as .csv
write.csv(Maynard_logCPM_dataset, file = here("Data", "processed_Data", 
                                               "Maynard_logCPM_dataset.csv"))
write.csv(Maynard_CPM_dataset, file = here("Data", "processed_Data", 
                                           "Maynard_CPM_dataset.csv"))

