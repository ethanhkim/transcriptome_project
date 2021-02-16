## Create transposed Maynard datasets

library(dplyr)
library(tidyverse)
library(magrittr)
library(conflicted)
library(vegan)

# Set conflicts
conflict_prefer('intersect', 'dplyr')

# Create common genelist between Maynard and He
He_genelist <- unique(He_DS1_averaged_by_layer$gene_symbol)
Maynard_genelist <- unique(Maynard_logCPM_averaged$gene_symbol)
common_genelist <- intersect(He_genelist, Maynard_genelist)


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

# Use non-averaged He and Maynard data

Maynard_transposed <- transpose_df(Maynard_logCPM_averaged, common_genelist)
He_transposed <- transpose_df(He_DS1_averaged_by_layer, common_genelist)

# Run Mantel tests

mantel_test <- function(df_1, df_2, no_of_permutations = 999) {
  
  df1_cor <- cor(df_1, df_1, method = "pearson", 
                 use = "complete.obs")
  df2_cor <- cor(df_2, df_2, method = "pearson", 
                 use = "complete.obs")
  
  mantel_test_result <- mantel(df1_cor, df2_cor,
                               method = "spearman",
                               permutations = no_of_permutations,
                               na.rm = T)
  
  return(mantel_test_result)
}

He_Maynard_mantel <- mantel_test(He_transposed, Maynard_transposed, 1)


testlist <- 

df2 %<>%
  filter(gene_symbol %in% df1_df2_genelist) %<>%
  arrange(gene_symbol) %<>%
  column_to_rownames(var = "gene_symbol") %<>%
  t()

df1_df2_cor <- cor(df1, df2, method = "pearson") 

df1_df2_cor <- df1_df2_cor %>%
  lower.tri(diag = TRUE) %>%
  as_tibble()
df1_test <- df[upper.tri(df1_df2_cor)]

head(df1_df2_cor)

  