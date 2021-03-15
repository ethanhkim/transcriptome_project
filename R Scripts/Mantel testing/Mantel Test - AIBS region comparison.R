## Mantel tests for between scRNA-seq data ##

# Load required libraries
library(dplyr)
library(tidyr)
library(data.table)
library(here)
library(magrittr)
library(edgeR)
library(conflicted)
library(tibble)
library(vegan)
library(parallel)
library(WGCNA)

# Declare conflicts
conflict_prefer("filter", "dplyr")

# Load AIBS_logCPM data
load(here("Data", "processed_data", "AIBS_logCPM.RData"))

# Separate by region and convert long to wide
separate_by_region <- function(all_regions_df, region) {
  
  df <- all_regions_df
  
  region_df <- df %>%
    filter(region_label == region) %>%
    spread(cortical_layer_label, log_mean_expression) %>%
    unite(metadata, c("gene_symbol", "region_label",
                      "class_label")) %>%
    column_to_rownames(var = "metadata")
  
  return(region_df)
}

# Create list of tibbles of all regions in layer(col) by gene(row) format
all_regions <- list()
for (region in unique(AIBS_logCPM$region_label)) {
  all_regions[[region]] <- separate_by_region(AIBS_logCPM, region)
}

# Function for mantel testing
mantel_test <- function(df_1, df_2, no_of_permutations = 999) {
  
  df1 <- df_1 %>% t()
  df2 <- df_2 %>% t()
  
  df1_cor <- WGCNA::cor(df_1, df_1, method = "pearson", 
                 use = "all.obs", nThreads = 12)
  df2_cor <- cor(df_2, df_2, method = "pearson", 
                 use = "all.obs", nThreads = 12)
  
  mantel_test_result <- mantel(df1_cor, df2_cor,
                               method = "spearman",
                               permutations = no_of_permutations,
                               na.rm = T,
                               parallel = 12)
  
  return(mantel_test_result)
}


# L1, L2, L3, L4, L5, L6a, L6b
S1_regions <- colnames(all_regions$S1lm)
# L1, L2, L3, L5, L6
M1_regions <- colnames(all_regions$M1lm)
# L1, L2, L3, L4, L5, L6a, L6b, WM
A1C_regions <- colnames(all_regions$A1C)
# L1, L2, L3, L4ab, L4c, L5, L6a, L6b
V1C_regions <- colnames(all_regions$V1C)
# L1, L2, L3, L4, L5, L6
MTG_regions <- colnames(all_regions$MTG)
# L1, L2, L5a, L5b, L6
CgG_regions <- colnames(all_regions$CgG)

