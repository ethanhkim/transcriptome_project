library(biomaRt)
library(EnsDb.Hsapiens.v79)
library(org.Hs.eg.db)
library(dplyr)
library(magrittr)
library(conflicted) # Easily manage conflicting libraries
library(here)

# Set conflicts
conflict_prefer('intersect', 'dplyr')
conflict_prefer("filter", "dplyr")
conflict_prefer("slice", "dplyr")

# Function to get biotypes
get_biotypes <- function(gene_symbol_col, gene_annotation_type) {

  metadata <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                    values = gene_symbol_col, mart = mart)
  if (gene_annotation_type == "ensembl_gene_id") {
    genes <- tibble(ensembl_gene_id = gene_symbol_col)
    metadata_filtered <- metadata %>%
      distinct(ensembl_gene_id, .keep_all = T)
    metadata_joined <- inner_join(metadata_filtered, genes, by = "ensembl_gene_id")
  } else {
    genes <- tibble(hgnc_symbol = gene_symbol_col)
    metadata_filtered <- metadata %>%
      distinct(hgnc_symbol, .keep_all = T)
    metadata_joined <- inner_join(metadata_filtered, genes, by = "hgnc_symbol")
  }
  
  metadata_joined %<>%
    group_by(gene_biotype) %>%
    summarize(n()) %>%
    arrange(desc(`n()`))
  
  return(metadata_joined)

}

# Load data
load(here("Data", "processed_data", "He_DS1_logCPM_dataset.Rdata"))
load(here("Data", "processed_data", "Maynard_logCPM_dataset.Rdata"))
load(here("Data", "processed_data", "He_DS1_logCPM_filtered_dataset.Rdata"))
load(here("Data", "processed_data", "Maynard_logCPM_filtered_dataset.Rdata"))


# CPM normalize, filter out CPM < 0.1
He_DS1_below_CPM <- He_DS1_sum_layer %>%
  # Add one to counts to avoid taking cpm of 0
  mutate_at(c("L1", "L2", "L3", "L4", "L5", "L6", "WM"), ~. +1) %>%
  column_to_rownames(var = "gene_symbol") %>%
  cpm() %>%
  as.data.frame() %>%
  # Retain gene symbols, as filter removes rownames
  rownames_to_column(var = "gene_symbol") %>%
  # Filter out samples of CPM < 0.1
  filter_at(vars(-gene_symbol), all_vars(. < .1))

# CPM normalize, filter out CPM < 0.1
He_DS1_above_CPM <- He_DS1_sum_layer %>%
  # Add one to counts to avoid taking cpm of 0
  mutate_at(c("L1", "L2", "L3", "L4", "L5", "L6", "WM"), ~. +1) %>%
  column_to_rownames(var = "gene_symbol") %>%
  cpm() %>%
  as.data.frame() %>%
  # Retain gene symbols, as filter removes rownames
  rownames_to_column(var = "gene_symbol") %>%
  # Filter out samples of CPM < 0.1
  filter_at(vars(-gene_symbol), all_vars(. > .1))

# Intersected gene lists
He_Maynard_logCPM_genes <- intersect(He_DS1_logCPM_dataset$gene_symbol, Maynard_logCPM_dataset$gene_symbol)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#Get gene biotypes

He_CPM_genes <- get_biotypes(He_count_matrix$gene_symbol, "hgnc_symbol")
He_low_CPM_genes <- get_biotypes(He_DS1_below_CPM$gene_symbol, "hgnc_symbol")
He_high_CPM_genes <- get_biotypes(He_DS1_above_CPM$gene_symbol, "hgnc_symbol")
Maynard_CPM_genes <- get_biotypes(Maynard_ensembl_list, "ensembl_gene_id")

He_Maynard_logCPM_biotypes <- get_biotypes(He_Maynard_logCPM_genes, "hgnc_symbol")








  
