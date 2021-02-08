## Create transposed Maynard datasets

library(dplyr)
library(tidyverse)
library(magrittr)
library(conflicted)
library(vegan)
library(here)
library(magrittr)
library(org.Hs.eg.db)

# Set conflicts
conflict_prefer('intersect', 'dplyr')

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

## DATA ##
# Use non-averaged He and Maynard data #

# Maynard data
Maynard_df <- read.csv(here('Data', 'raw_data', 'Maynard et al',
                            'layer_level_data.csv'),
                       stringsAsFactors = FALSE) %>%
  # Remove X's from col names
  rename_at(vars(contains('X')), funs(sub('X', '', .)))

# Map Ensembl to HGNC
Maynard_df$gene_symbol <- mapIds(org.Hs.eg.db, 
                                 keys = Maynard_df$Ensembl_ID, 
                                 keytype = "ENSEMBL", column="SYMBOL")
# Format Maynard data
Maynard_df %<>% 
  select(-"Ensembl_ID") %>% 
  select("gene_symbol", everything()) %>%
  filter(!is.na(gene_symbol)) %>%
  distinct(gene_symbol, .keep_all = TRUE)


# He data - run 1. Human - Preprocessing - He et al.R script
# use He_DS1_Human, removing S17 (white matter)
He_df <- He_DS1_Human %>%
  select(-S17)

# Create common genelist between Maynard and He
He_genelist <- unique(He_df$gene_symbol)
Maynard_genelist <- unique(Maynard_df$gene_symbol)
common_genelist <- intersect(He_genelist, Maynard_genelist)


# Transpose data to be in the shape of cuts by genes to correlate #
Maynard_transposed <- transpose_df(Maynard_df, common_genelist)
He_transposed <- transpose_df(He_df, common_genelist)

# Run Mantel tests

mantel_test <- function(df_1, df_2, no_of_permutations = 999) {
  
  df1_cor <- cor(df_1, df_1, method = "pearson")
  df2_cor <- cor(df_2, df_2, method = "pearson")
  
  df1_cor[is.na(df1_cor)] = 0
  df2_cor[is.na(df2_cor)] = 0
  
  mantel_test_result <- mantel(df1_cor, df2_cor,
                               method = "spearman",
                               permutations = no_of_permutations)
  
  return(mantel_test_result)
}

He_Maynard_mantel <- mantel_test(He_transposed, Maynard_transposed, 1)

# Testing Grounds #
Zeng_genelist <- Zeng_dataset_updated$gene_symbol

test_genelist <- intersect(common_genelist, Zeng_genelist)
test_genelist_markers <- c("NDNF", "CHRNA7", "CNR1",
                           'CXCL14', "DISC1", "INPP4B",
                           "RELN")

Maynard_test <- Maynard_df %>%
  select(-contains("WM"))

He_subset <- transpose_df(He_df, test_genelist) %>%
  replace(is.na(.), 0)
Maynard_subset <- transpose_df(Maynard_df, test_genelist_markers)


He_test_cor <- cor(He_subset, He_subset, method = "pearson")
Maynard_test_cor <- cor(Maynard_subset, Maynard_subset,
                        method = "pearson")

He_test_cor[is.na(He_test_cor)] <- 0
Maynard_test_cor[is.na(Maynard_test_cor)] <- 0


mantel_test(He_subset, Maynard_subset, no_of_permutations = 1)


test_cor <- cor(He_transposed, He_transposed, method = "pearson")
test_cor_2 <- cor(Maynard_transposed, Maynard_transposed, method = "pearson")