## Mantel test ##

library(tidyverse)
library(magrittr)
library(here)
library(data.table)
library(vegan) # required for mantel()
library(conflicted)

# Set conflicts

conflict_prefer('filter', 'dplyr')
conflict_prefer('select', 'dplyr')

# Import data #

# Maynard - logCPM data #
Maynard_logCPM_averaged <- fread(
  here("Data", "processed_data",
       "Maynard_logCPM_averaged.csv"))
Maynard_logCPM_averaged %<>% select(-V1)


# AIBS data #
AIBS_MTG <- fread("Data/Allen/MTG_df_01_21.csv") %>%
  # Remove col number column
  select(-V1)

# Separate by cell type and widen
widen_by_cell_type <- function(AIBS_df, cell_type, location = 'local') {
  df <- AIBS_df %>% 
    filter(class_label == cell_type)
  
  if (location == 'local') {
    df %<>% pivot_wider(
      names_from = cortical_layer_label, 
      values_from = mean_expression
    ) %>%
      select(-class_label)
  } else if (location == 'SCC') {
      df %<>% spread(
        key = cortical_layer_label, 
        value = mean_expression) %>%
      select(-class_label)
  }
  
  return(df)
}

# Filter out for common genes

AIBS_genelist <- unique(AIBS_MTG$gene_symbol)
Maynard_genelist <- unique(Maynard_logCPM_averaged$gene_symbol)
common_genelist <- intersect(AIBS_genelist, Maynard_genelist)


# Transpose to layer (row) by gene (column) for correlation
transpose_df <- function(subset_expression_df) {
  df <- subset_expression_df
  
  df %<>% column_to_rownames(var = "gene_symbol") %>%
    t() %>%
    as.matrix()
  
  return(df)
}

# Matrices to compare
AIBS_GABA <- widen_by_cell_type(AIBS_MTG, "GABAergic", 'SCC') %>%
  filter(gene_symbol %in% common_genelist) %>%
  transpose_df()
AIBS_GLUT <- widen_by_cell_type(AIBS_MTG, "Glutamatergic", "SCC") %>%
  filter(gene_symbol %in% common_genelist) %>%
  transpose_df()
AIBS_NONN <- widen_by_cell_type(AIBS_MTG, "Non-neuronal", "SCC") %>%
  filter(gene_symbol %in% common_genelist) %>%
  transpose_df()

Maynard_df <- Maynard_logCPM_averaged %>% 
  filter(gene_symbol %in% common_genelist) %>%
  arrange(gene_symbol) %>%
  transpose_df()

# Perform Mantel test

mantel_test <- function(df_1, df_2) {
  # Create correlation matrix for each dataset 
  df1_cor <- cor(df_1, df_1, use = "complete.obs")
  df2_cor <- cor(df_2, df_2, use = "complete.obs")
  
  mantel(df1_cor, df2_cor, method = "pearson", na.rm = T)
  
}

mantel_GABA <- mantel_test(AIBS_GABA, Maynard_df)
mantel_GLUT <- mantel_test(AIBS_GLUT, Maynard_df)
mantel_NONN <- mantel_test(AIBS_NONN, Maynard_df)

