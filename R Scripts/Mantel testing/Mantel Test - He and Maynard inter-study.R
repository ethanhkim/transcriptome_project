## Create transposed Maynard datasets ##

library(dplyr)
library(magrittr)
library(conflicted) # Easily manage conflicting libraries
library(vegan) # Package for mantel test
library(here)
library(parallel) # Parallelize Mantel and WGCNA corr()
library(WGCNA) # Faster corr() than base

# Set conflicts
conflict_prefer('intersect', 'dplyr')
conflict_prefer("filter", "dplyr")
conflict_prefer("slice", "dplyr")

# Function to perform Mantel testing
mantel_test <- function(He_df, Maynard_df, no_of_perm = 1) {
  # Create common gene list between Maynard and He
  create_common_genelist <- function(He_df, Maynard_df) {
    
    if ((isFALSE(duplicated(He_df))) == 
        (isFALSE(duplicated(Maynard_df)))) {
      
      common_genelist <- intersect(He_df$gene_symbol, 
                                   Maynard_df$gene_symbol)
    }
    return(common_genelist)
  }
  
  # Create transposed dataframe
  transpose_df <- function(df, genelist) {
    df %<>%
      as.data.frame() %>%
      filter(gene_symbol %in% genelist) %>%
      slice(match(genelist, gene_symbol)) %>%
      column_to_rownames(var = "gene_symbol") %>%
      t()
    return(df)
  }
  
  # Make common genelist
  common_genelist <- create_common_genelist(He_df,
                                            Maynard_df)
  
  print(length(common_genelist))
  
  # Transpose data
  He_transposed <- transpose_df(He_df, common_genelist)
  Maynard_transposed <- transpose_df(Maynard_df, common_genelist)
  
  # Create correlation matrices
  He_corr_matrix <- WGCNA::cor(He_transposed, method = "pearson",
                               use = "all.obs", nThreads = 12)
  Maynard_corr_matrix <- WGCNA::cor(Maynard_transposed, method = "pearson",
                                    use = "all.obs", nThreads = 12)
  
  mantel(He_corr_matrix, Maynard_corr_matrix, 
         permutations = no_of_perm,
         parallel = 12, na.rm = T)
}

# Load in data
load(here("Data", "processed_data", "He_DS1_logCPM_dataset.Rdata"))
load(here("Data", "processed_data", "Maynard_logCPM_dataset.Rdata"))
load(here("Data", "processed_data", "He_DS1_logCPM_filtered_dataset.Rdata"))
load(here("Data", "processed_data", "Maynard_logCPM_filtered_dataset.Rdata"))
load(here("Data", "processed_data", "Zeng_dataset_long.Rdata"))
Zeng_marker_genes <- unique(Zeng_dataset_long$gene_symbol)


## Mantel tests b/w He and Maynard ----

# Non-logCPM data
# Mantel r: 0.4079
mantel_test(He_DS1_averaged_by_layer, Maynard_dataset_average)

# He logCPM vs Maynard non-logCPM ----
# Mantel r: 0.4392
mantel_test(He_DS1_logCPM_dataset, Maynard_dataset_average)

# He non-logCPM vs Maynard logCPM ----
# Mantel r: 0.4646
mantel_test(He_DS1_averaged_by_layer, Maynard_logCPM_dataset)

# He logCPM vs Maynard logCPM ----
# Mantel r: 0.5166, p < 0.001
mantel_test(He_DS1_logCPM_dataset, Maynard_logCPM_dataset)

# He CPM vs Maynard CPM - does the score get better or worse?
# Mantel r: 0.4954; marginally worse but not by much
mantel_test(He_DS1_CPM_dataset, Maynard_CPM_dataset)

# He logCPM vs Maynard logCPM filtered for CPM > 0.1
# Mantel r: 0.5731
mantel_test(He_DS1_logCPM_filtered_dataset, Maynard_logCPM_filtered_dataset)

# Subset He and Maynard by Zeng markers
He_Zeng_subset <- He_DS1_logCPM_dataset %>%
  filter(gene_symbol %in% Zeng_marker_genes)
Maynard_Zeng_subset <- Maynard_logCPM_dataset %>%
  filter(gene_symbol %in% Zeng_marker_genes)
# He logCPM vs Maynard logCPM - subset by Zeng
# Mantel r: 0.683, p < 0.001
mantel_test(He_Zeng_subset, Maynard_Zeng_subset)

# Subset He and Maynard by Zeng markers, with CPM > 0.1
He_Zeng_subset_filtered <- He_DS1_logCPM_filtered_dataset %>%
  filter(gene_symbol %in% Zeng_marker_genes)
Maynard_Zeng_subset_filtered <- Maynard_logCPM_filtered_dataset %>%
  filter(gene_symbol %in% Zeng_marker_genes)
# He logCPM vs Maynard logCPM - subset by Zeng
# Mantel r: 0.683, p < 0.001
mantel_test(He_Zeng_subset_filtered, Maynard_Zeng_subset_filtered)





