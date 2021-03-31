## Allen vs. Maynard vs. He ##

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

# Load in data
load(here("Data", "processed_data", "He_DS1_logCPM_dataset.Rdata"))
load(here("Data", "processed_data", "Maynard_logCPM_dataset.Rdata"))
load(here("Data", "processed_data", "He_DS1_logCPM_filtered_dataset.Rdata"))
load(here("Data", "processed_data", "Maynard_logCPM_filtered_dataset.Rdata"))
load(here("Data", "processed_data", "Allen_logCPM_dataset.Rdata"))
load(here("Data", "processed_data", "Allen_logCPM_filtered_dataset.Rdata"))

# Zeng markers
load(here("Data", "processed_data", "Zeng_dataset_long.Rdata"))
Zeng_marker_genes <- unique(Zeng_dataset_long$gene_symbol)

# Cell-type specific AIBS
Allen_GABA <- Allen_logCPM_dataset %>%
  filter(class_label == "GABAergic") %>%
  select(-class_label)
Allen_GLUT <- Allen_logCPM_dataset %>%
  filter(class_label == "Glutamatergic") %>%
  select(-class_label)
Allen_NONN <- Allen_logCPM_dataset %>%
  filter(class_label == "Non-neuronal") %>%
  select(-class_label)

## Mantel tests b/w He, Maynard, Allen ##

# He logCPM vs Maynard logCPM
# Mantel r: 0.5166, p < 0.001, n = 18,206
mantel_test(He_DS1_logCPM_dataset, Maynard_logCPM_dataset)

# He logCPM vs Allen logCPM - GABAergic
# Mantel r: 0.5166, p < 0.001, n = 30,744
mantel_test(He_DS1_logCPM_dataset, Allen_GABA)

# He logCPM vs Allen logCPM - Glutamatergic
# Mantel r: 0.5166, p < 0.001, n = 30,744
mantel_test(He_DS1_logCPM_dataset, Allen_GLUT)

# He logCPM vs Allen logCPM - Non-neuronal
# Mantel r: 0.5166, p < 0.001, n = 30,744
mantel_test(He_DS1_logCPM_dataset, Allen_NONN)


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
