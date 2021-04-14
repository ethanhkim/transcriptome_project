## Mantel Test - Layer-to-layer comparisons ##

# Load libraries
library(dplyr)
library(magrittr)
library(tibble)
library(conflicted) # Easily manage conflicting libraries
library(vegan) # Package for mantel test
library(here)
library(parallel) # Parallelize Mantel and WGCNA corr()
library(WGCNA) # Faster corr() than base
library(gplots)

# Set conflicts
conflict_prefer('intersect', 'dplyr')
conflict_prefer("filter", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("cor", "WGCNA")

# Source mantel test function
source(here("R Scripts", "Mantel testing", "mantel_test.R"))

# Load in data
load(here("Data", "processed_data", "He_DS1_logCPM_dataset.Rdata"))
load(here("Data", "processed_data", "He_DS1_logCPM_filtered_dataset.Rdata"))
load(here("Data", "processed_data", "Maynard_logCPM_dataset.Rdata"))
load(here("Data", "processed_data", "Maynard_logCPM_filtered_dataset.Rdata"))
load(here("Data", "processed_data", "Allen_gene_logCPM_dataset.Rdata"))
load(here("Data", "processed_data", "Allen_gene_logCPM_filtered_dataset.Rdata"))

# Format data
Allen_gene_logCPM_dataset %<>% column_to_rownames(var = "gene_symbol") 
Allen_gene_logCPM_filtered_dataset %<>% column_to_rownames(var = "gene_symbol")
He_DS1_logCPM_dataset %<>% column_to_rownames(var = "gene_symbol")
He_DS1_logCPM_filtered_dataset %<>% column_to_rownames(var = "gene_symbol")
Maynard_logCPM_dataset %<>% column_to_rownames(var = "gene_symbol")
Maynard_logCPM_filtered_dataset %<>% column_to_rownames(var = "gene_symbol")


# Create heatmaps for each of the layer-to-layer correlations
layer_cor_heatmaps <- function(source_data_1, source_data_2) {
  cor_df <- cor(source_data_1, source_data_2,
                use = 'pairwise.complete.obs')
  heatmap.2(cor_df)
}

layer_cor_heatmaps(He_DS1_logCPM_dataset, He_DS1_logCPM_dataset)
layer_cor_heatmaps(He_DS1_logCPM_dataset, Maynard_logCPM_dataset)
layer_cor_heatmaps(He_DS1_logCPM_filtered_dataset)

layer_cor_heatmaps(Maynard_logCPM_dataset)
layer_cor_heatmaps(Maynard_logCPM_filtered_dataset)

layer_cor_heatmaps(Allen_gene_logCPM_dataset)
layer_cor_heatmaps(Allen_gene_logCPM_filtered_dataset)

