## Create transposed Maynard datasets

library(dplyr)
library(tidyverse)
library(magrittr)


#Function to create transposed dataset

transpose_Maynard_dataset <- function(x) {
  as.data.frame(x) %>%
    filter(gene_symbol %in% Maynard_common_gene_list) %>%
    slice(match(common_gene_list, gene_symbol)) %>%
    column_to_rownames(var = "gene_symbol") %>%
    t()
}

Maynard_dataset_average_transposed <- transpose_Maynard_dataset(Maynard_dataset_average)
Maynard_dataset_transposed <- transpose_Maynard_dataset(Maynard_dataset)
  
# Create correlation matrices

He_values_matrix <- cor(He_values_transposed, He_values_transposed, method = "pearson")
Maynard_values_matrix <- cor(Maynard_dataset_transposed, Maynard_dataset_transposed, method = "pearson")

results <- mantel(He_values_matrix, Maynard_values_matrix, method = "spearman", permutations = 999) %>%
  print()

Mantel_test_result_statistic <- results$statistic


df1 <- Maynard_dataset_average 
df2 <- He_DS1_averaged_by_layer

df1_df2_genelist <- intersect(df1$gene_symbol, df2$gene_symbol)

df1 %<>%
  filter(gene_symbol %in% df1_df2_genelist) %<>%
  arrange(gene_symbol) %<>%
  column_to_rownames(var = "gene_symbol") %<>%
  t()

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


  