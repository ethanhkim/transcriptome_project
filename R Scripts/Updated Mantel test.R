## Mantel test ##


library(tidyverse)
library(magrittr)
library(here)
library(data.table)
library(vegan) # required for mantel()


# Import data #

# Maynard - logCPM data #

Maynard_df <- Maynard_logCPM_averaged

# AIBS data #
AIBS_MTG <- fread("Data/Allen/MTG_df_01_21.csv") %>%
  # Remove col number column
  select(-V1)

# Separate by cell type and widen
widen_by_cell_type <- function(AIBS_df, cell_type) {
  df <- AIBS_df
  
  df %<>% filter(class_label == cell_type) %>%
    pivot_wider(
      names_from = cortical_layer_label, values_from = mean_expression
    ) %>%
    select(-class_label)
  
  return(df)
}

# Filter out for common genes

AIBS_genelist <- unique(AIBS_MTG$gene_symbol)
Maynard_genelist <- unique(Maynard_logCPM_averaged$gene_symbol)

common_genelist <- intersect(AIBS_genelist, Maynard_genelist)

Maynard_df %<>% filter(gene_symbol %in% common_genelist)
AIBS_MTG %<>% filter(gene_symbol %in% common_genelist)

# Transpose to layer (row) by gene (column) for correlation
transpose_df <- function(subset_expression_df) {
  df <- subset_expression_df
  
  df %<>% column_to_rownames(var = "gene_symbol") %>%
    t() %>%
    as.matrix()
  
  return(df)
}

AIBS_GABA <- widen_by_cell_type(AIBS_MTG, "GABAergic") %>%
  transpose_df()
Maynard_df %<>% transpose_df()

# Perform Mantel test

mantel_test <- function(df_1, df_2) {
  
  df1_cor <- cor(df_1, df_1, use = "complete.obs")
  df2_cor <- cor(df_2, df_2, use = "complete.obs")
  
  mantel(df1_cor, df2_cor, method = "pearson", na.rm = T)
  
}

test <- mantel_test(AIBS_GABA, Maynard_df)

test1 <- cor(AIBS_GABA, AIBS_GABA, use = "complete.obs")
test2 <- cor(Maynard_df, Maynard_df)
mantel(test1, test2)
