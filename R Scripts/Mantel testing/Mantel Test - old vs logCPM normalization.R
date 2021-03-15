## Mantel tests between old bulk-tissue vs. logCPM normalized bulk-tissue ##

library(dplyr)
library(tibble)
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
mantel_test <- function(df_1, df_2, no_of_perm = 1) {
  # Create common gene list between Maynard and He
  create_common_genelist <- function(df_1, df_2) {
    
    if ((isFALSE(duplicated(df_1))) == 
        (isFALSE(duplicated(df_2)))) {
      
      common_genelist <- intersect(df_1$gene_symbol, 
                                   df_2$gene_symbol)
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
  common_genelist <- create_common_genelist(df_1,
                                            df_2)
  
  # Transpose data
  He_transposed <- transpose_df(df_1, common_genelist)
  Maynard_transposed <- transpose_df(df_2, common_genelist)
  
  # Create correlation matrices
  He_corr_matrix <- WGCNA::cor(He_transposed, method = "pearson",
                               use = "all.obs", nThreads = 12)
  Maynard_corr_matrix <- WGCNA::cor(Maynard_transposed, method = "pearson",
                                    use = "all.obs", nThreads = 12)
  
  mantel(He_corr_matrix, Maynard_corr_matrix, 
         permutations = no_of_perm,
         parallel = 12, na.rm = T)
}

# Load in data #

# Old bulk-tissue
load(here("Data", "processed_data", "He_DS1_averaged_by_layer.Rdata"))
load(here("Data", "processed_data", "Maynard_dataset_average.Rdata"))
# logCPM bulk-tissue
load(here("Data", "processed_data", "He_DS1_logCPM_dataset.Rdata"))
load(here("Data", "processed_data", "Maynard_logCPM_dataset.Rdata"))

# old He average vs. logCPM
# Mantel r: 0.6446
mantel_test(He_DS1_averaged_by_layer, He_DS1_logCPM_dataset)
# old Maynard average vs. logCPM
# Mantel r: 0.8379
mantel_test(Maynard_dataset_average, Maynard_logCPM_dataset)


# Conclusion: looks like the old and new datasets are quite similar.
# Should be safe to use this data in the webapp.