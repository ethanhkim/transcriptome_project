# Load necessary libraries #


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("LieberInstitute/spatialLIBD")
library(spatialLIBD)
library(tidyverse)
library(org.Hs.eg.db)
library(magrittr)
library(data.table)
library(edgeR)
library(biomaRt)
library(stringr)
library(conflicted)
library(here)


# Set conflicts #
conflict_prefer('select', 'dplyr')
conflict_prefer('filter', 'dplyr')
conflict_prefer('cpm', 'edgeR')
conflict_prefer('rename', 'dplyr')

# Function to create columns of mean values across individuals for Maynard
create_mean_col <- function(Maynard_normalized_subset, layer_label) {
  df <- Maynard_normalized_subset
  layer <- layer_label
  
  df %<>% select(contains(layer)) %>%
    mutate(mean = rowMeans(.)) %>%
    pull(mean)
  
  return(df)
}

## Maynard normalization ##

# Fetch data from spatialLIBD
sce_layer <- fetch_data(type = "sce_layer")
# Get raw UMI counts per gene
Maynard_dataset <- sce_layer@assays@data@listData$counts %>%
  as_tibble(rownames=NA)
# Get gene Ensembl ID
Maynard_ensembl_list <- sce_layer@rowRanges@ranges@NAMES

# Normalize Maynard UMI counts with log2CPM
Maynard_logCPM <- cpm(Maynard_dataset, log = T) %>%
  as.data.frame() %>%
  rownames_to_column(var='Ensembl_ID')

# Convert Ensembl to HGNC, filter out non-viable HGNC
Maynard_logCPM$gene_symbol <- mapIds(org.Hs.eg.db, keys = Maynard_logCPM$Ensembl_ID, keytype = "ENSEMBL", column="SYMBOL")
Maynard_logCPM %<>%
  dplyr::select(-"Ensembl_ID") %>% 
  dplyr::select("gene_symbol", everything()) %>%
  filter(!is.na(gene_symbol)) %>%
  distinct(gene_symbol, .keep_all = TRUE)

# Select individuals with all layers
Maynard_logCPM_subset <- Maynard_logCPM %>%
  select("gene_symbol", contains('151507'), contains('151508'), contains('151509'),
         contains('151673'), contains('151674'), contains('151675'), contains('151676')) %>%
  column_to_rownames(var = "gene_symbol")

# Create list of mean columns
Maynard_mean_col_list <- list()
for (label in c("Layer1", "Layer2", "Layer3", "Layer4", "Layer5", "Layer6")) {
  Maynard_mean_col_list[[label]] <- create_mean_col(Maynard_logCPM_subset, label)
}
# Bind list and convert to dataframe
Maynard_logCPM_averaged <- as.data.frame(do.call(rbind, Maynard_mean_col_list)) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var = 'gene_symbol')

# Write Maynard normalized data

save(Maynard_logCPM_averaged, file = here("Data", "processed_data", 
                                          "Maynard_logCPM_averaged.Rdata"))
write.csv(Maynard_logCPM_averaged, file = here("Data", "processed_Data", 
                                               "Maynard_logCPM_averaged.csv"))


## He et al Normalization ##

# http://seqanswers.com/forums/showthread.php?t=59202 
# Converting RPKM values to logCPM values:
# RPKM = 2^(logCPM - log2(genelength))
# log2(RPKM) = logCPM - log2(genelength)
# log2(RPKM) + log2(genelength) = logCPM
# log2(RPKM * genelength) = logCPM

## DS1 ##
He_DS1_Human <- read.table(here("Data", "raw_data", "He et al", "tab.average_RPKM_section_human.tsv"), header = T)
# Truncate .xx from Ensembl ID's
He_DS1_Human$ENSEMBL_ID <- gsub("\\.\\d+$", "", He_DS1_Human$ENSEMBL_ID)
# Get gene length of the Ensembl ID's
He_Ensembl_ID_genelength <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_coords <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "start_position", "end_position"),
                     mart = He_Ensembl_ID_genelength)
gene_coords$size <- gene_coords$end_position - gene_coords$start_position
gene_coords %<>% filter(ensembl_gene_id %in% He_DS1_Human$ENSEMBL_ID) %>%
  rename(ENSEMBL_ID = ensembl_gene_id, gene_symbol = hgnc_symbol) %>%
  select(ENSEMBL_ID, gene_symbol, size)

# Set 0's to -Inf
He_DS1_Human[He_DS1_Human == 0] <- -Inf

He_DS1_Human %<>%
  inner_join(gene_coords, by = 'ENSEMBL_ID') %>%
  select(gene_symbol, everything()) %>%
  select(-ENSEMBL_ID) %>%
  filter(!is.na(gene_symbol)) %>%
  distinct(gene_symbol, .keep_all = TRUE)

gene_length_vector <- He_DS1_Human$size

He_DS1_Human %<>%
  select(-size) %>%
  column_to_rownames(var = 'gene_symbol')

He_DS1_Human <- He_DS1_Human * gene_length_vector
He_DS1_log2 <- log(He_DS1_Human, 2) 
He_DS1_log_scale <- scale(He_DS1_log)


# Multiply counts by gene length

# Log multiplied counts to get CPM



## Getting AIBS normalized data ##

# Function to separate out by cell type

separate_by_cell_type <- function(AIBS_MTG_df, cell_type) {
  df <- AIBS_MTG_df

  df %<>% 
    # Filter by specified cell type
    filter(class_label == cell_type) %>%
    # Widen to gene by layer
    pivot_wider(names_from = cortical_layer_label,
                values_from = mean_expression) %>%
    # Gene symbols to row name
    column_to_rownames(var = 'gene_symbol')
  
  return(df)
}

<<<<<<< HEAD
# Read in processed csv
=======
# Read in processed csv, produced from scRNA-seq CPM script
>>>>>>> ec361f9ef70abd13a78bc178bb3f36951c54a147

AIBS_MTG <- fread(here('Data', 'Allen', 'MTG_df_01_21.csv')) %>%
  # remove col number column
  select(-V1)

AIBS_MTG_cell_type <- list()
for (i in c("GABAergic", "Glutamatergic", "Non-neuronal")) {
  AIBS_MTG_cell_type[[i]] <- separate_by_cell_type(AIBS_MTG, i)
}

AIBS_MTG_all_cell_types <- rbind(AIBS_MTG_cell_type$GABAergic, 
                               AIBS_MTG_cell_type$Glutamatergic, 
                               AIBS_MTG_cell_type$`Non-neuronal`)

# Function to create cell
AIBS_MTG_GABA <- AIBS_MTG_all_cell_types %>%
  filter(class_label == 'GABAergic') %>%
  select(-class_label) %>%
  t() %>%
  scale() %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var = 'gene_symbol') %>%
  pivot_longer(cols = L1:L6,
               names_to = 'cortical_layer_label',
               values_to = "mean_expression") %>%
  add_column(class_label = 'GABAergic')



