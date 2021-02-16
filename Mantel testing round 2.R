# Source data:
Maynard_df <- read.csv(here('Data', 'raw_data', 'Maynard et al',
                                 'layer_level_data.csv'),
                            stringsAsFactors = FALSE) %>%
  # Remove X's from col names
  rename_at(vars(contains('X')), funs(sub('X', '', .)))

He_df <- He_DS1_Human

# Create transposed dataframe
transpose_df <- function(df, genelist) {
  df %<>%
    as.data.frame() %>%
    filter(gene_symbol %in% genelist) %>%
    slice(match(genelist, gene_symbol)) %>%
    column_to_rownames(var = "gene_symbol") %>%
    t() %>%
    replace(is.na(.), 0)
  return(df)
}

# Genelist to use:
test_genelist <- intersect(common_genelist, Zeng_genelist)
test_genelist_markers <- c("NDNF", "CHRNA7", "CNR1",
                           'CXCL14', "DISC1", "INPP4B",
                           "RELN")
common_genelist <- intersect(He_df$gene_symbol, Maynard_df$gene_symbol)

## Subset 1: No modifications ##
# Maynard - not quite correct.
Maynard_subset_1 <- transpose_df(Maynard_dataset, test_genelist_markers)
Maynard_cor_1 <- cor(Maynard_subset_1, Maynard_subset_1,
                        method = "pearson") %>%
  as_tibble(rownames = NA)

# H - this is the correct dataframe used to replicate earliest figure #
He_subset_1 <- He_df %>% select(-S17) %>%
  transpose_df(genelist = test_genelist_markers)
He_cor_1 <- cor(He_subset_1, He_subset_1,
                method = "pearson") %>%
  as_tibble(rownames = NA)


## Subset 2 - Maynard: ##
# Maynard - this is the dataframe used for earliest figure
Maynard_subset <- Maynard_dataset %>%
  select(-contains("WM")) %>%
  transpose_df(genelist = test_genelist_markers)
Maynard_cor <- cor(Maynard_subset, Maynard_subset,
                     method = "pearson") %>%
  as_tibble(rownames = NA)

# DF's without WM, selected for Zeng et al genes
He_subset <- He_df %>% select(-S17) %>%
  transpose_df(genelist = test_genelist)
Maynard_subset <- Maynard_df %>% 
  select(-contains("WM")) %>%
  transpose_df(genelist = test_genelist)

He_cor <- cor(He_subset, He_subset, method = "pearson")
He_cor[is.na(He_cor)] <- 0
Maynard_cor <- cor(Maynard_subset, Maynard_subset, method = "pearson")

mantel(He_cor, Maynard_cor, method = "spearman", permutations = 1)

# DF's with WM, selected for Zeng et al genes
He_subset <- He_df %>%
  transpose_df(genelist = test_genelist)
Maynard_subset <- Maynard_df %>% 
  transpose_df(genelist = test_genelist)

He_cor <- cor(He_subset, He_subset, method = "pearson")
He_cor[is.na(He_cor)] <- 0
Maynard_cor <- cor(Maynard_subset, Maynard_subset, method = "pearson")

mantel(He_cor, Maynard_cor, method = "spearman", permutations = 1)

# DF's with WM, selected for common genes across Maynard and He
He_subset <- He_df %>%
  transpose_df(genelist = common_genelist)
Maynard_subset <- Maynard_df %>% 
  transpose_df(genelist = common_genelist)

He_cor <- cor(He_subset, He_subset, method = "pearson")
He_cor[is.na(He_cor)] <- 0
Maynard_cor <- cor(Maynard_subset, Maynard_subset, method = "pearson")

mantel(He_cor, Maynard_cor, method = "spearman", permutations = 1)

