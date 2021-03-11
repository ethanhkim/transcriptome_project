## Create transposed Maynard datasets

library(dplyr)
library(magrittr)
library(conflicted)
library(vegan)
library(here)
library(vegan)
library(parallel)
library(WGCNA)

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

# Load in data
load(here("Data", "processed_data", "He_DS1_logCPM_dataset.Rdata"))
load(here("Data", "processed_data", "Maynard_logCPM_dataset.Rdata"))


## Mantel tests ##

# Old non-logCPM averaged data for He and Maynard
mantel_test(He_DS1_averaged_by_layer, Maynard_dataset_average)
# He logCPM vs old Maynard
mantel_test(He_DS1_logCPM_dataset, Maynard_dataset_average)
# old He vs Maynard logCPM
mantel_test(He_DS1_averaged_by_layer, Maynard_logCPM_dataset)
# He logCPM vs Maynard logCPM
mantel_test(He_DS1_logCPM_dataset, Maynard_logCPM_dataset, no_of_perm = 1000)


  