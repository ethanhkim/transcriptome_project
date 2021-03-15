## Mantel tests between old bulk-tissue vs. logCPM normalized bulk-tissue ##

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

# Load in data #

# Old bulk-tissue
load(here("Data", "processed_data", "He_DS1_averaged_by_layer.Rdata"))
load(here("Data", "processed_data", "Maynard_dataset_average.Rdata"))
# logCPM bulk-tissue
load(here("Data", "processed_data", "He_DS1_logCPM_dataset.Rdata"))
load(here("Data", "processed_data", "Maynard_logCPM_dataset.Rdata"))


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