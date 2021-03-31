#### Mantel Test function used for all Mantel Test scripts ##


# Load required libraries

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
  common_genelist <- create_common_genelist(df_1, df_2)
  
  print(length(common_genelist))
  
  # Transpose data
  df_1_transposed <- transpose_df(df_1, common_genelist)
  df_2_transposed <- transpose_df(df_2, common_genelist)
  
  # Create correlation matrices
  df_1_corr_matrix <- WGCNA::cor(df_1_transposed, method = "pearson",
                                 use = "pairwise.complete.obs", nThreads = 12)
  df_2_corr_matrix <- WGCNA::cor(df_2_transposed, method = "pearson",
                                 use = "pairwise.complete.obs", nThreads = 12)
  
  mantel(df_1_corr_matrix, df_2_corr_matrix, 
         permutations = no_of_perm,
         parallel = 12, na.rm = T)
}

