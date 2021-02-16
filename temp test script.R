# Using preprocessCore or limma

library(preprocessCore)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(conflicted)

conflict_prefer('filter', 'dplyr')

## Functions ##

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

## Data and lists ##

# He data - stock: run 1 - He et al Preprocessing

# List of genes for He:
He_genes <- He_DS1_Human %>%
  rownames_to_column(var = "gene_symbol") %>%
  pull(var = "gene_symbol")

# He quantile-normalized
# Uses He_DS1_Human after running 1 - He et al Preprocessing
He_matrix <- He_DS1_Human %>%
  as.matrix()

He_ranked_matrix <- normalize.quantiles(He_matrix)

# Maynard data - stock
Maynard_df <- read.csv(here('Data', 'raw_data', 'Maynard et al',
                            'layer_level_data.csv'),
                       stringsAsFactors = FALSE) %>%
  # Remove X's from col names
  rename_at(vars(contains('X')), funs(sub('X', '', .)))

# Maynard data - logCPM
Maynard_logCPM <- Maynard_df %>%
  column_to_rownames(var = "gene_symbol") %>%
  cpm(log = T) %>%
  as.data.frame() %>%
  rownames_to_column(var='gene_symbol')

# Genelist to use:
common_genelist <- intersect(He_df$gene_symbol, Maynard_df$gene_symbol)
test_genelist <- intersect(common_genelist, Zeng_genelist)
test_genelist_markers <- c("NDNF", "CHRNA7", "CNR1",
                           'CXCL14', "DISC1", "INPP4B",
                           "RELN")



# Testing Mantel correlations # 
# DF's with WM, selected for Zeng et al genes

# Set 1: He (stock) + Maynard (stock)

He_subset <- He_DS1_Human %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_symbol") %>%
  transpose_df(genelist = test_genelist)
Maynard_subset <- Maynard_df %>%
  transpose_df(genelist = test_genelist)

He_cor <- cor(He_subset, He_subset, method = "pearson")
He_cor[is.na(He_cor)] <- 0
Maynard_cor <- cor(Maynard_subset, Maynard_subset, method = "pearson")

mantel(He_cor, Maynard_cor, method = "spearman", permutations = 1)
# Mantel statistic r: 0.6257

# Set 2: He (quantile) + Maynard (stock)

He_subset <- He_ranked_matrix %>%
  as.data.frame() %>% add_column(gene_symbol = He_genes) %>%
  transpose_df(genelist = test_genelist)
Maynard_subset <- Maynard_df %>% 
  transpose_df(genelist = test_genelist)

He_cor <- cor(He_subset, He_subset, method = "pearson")
He_cor[is.na(He_cor)] <- 0
Maynard_cor <- cor(Maynard_subset, Maynard_subset, method = "pearson")

mantel(He_cor, Maynard_cor, method = "spearman", permutations = 1)
# Mantel statistic r: 0.6151

# Set 2: He (quantile) + Maynard (logCPM)

He_subset <- He_ranked_matrix %>%
  as.data.frame() %>% add_column(gene_symbol = He_genes) %>%
  transpose_df(genelist = test_genelist)
Maynard_subset <- Maynard_logCPM %>% 
  transpose_df(genelist = test_genelist)

He_cor <- cor(He_subset, He_subset, method = "pearson")
He_cor[is.na(He_cor)] <- 0
Maynard_cor <- cor(Maynard_subset, Maynard_subset, method = "pearson")

mantel(He_cor, Maynard_cor, method = "spearman", permutations = 1)
# Mantel statistic r: 0.5882

# Set 3: He (quantile + Z-score) + Maynard (logCPM + Z-score)

He_subset <- He_ranked_matrix %>%
  as.data.frame() %>% add_column(gene_symbol = He_genes) %>%
  column_to_rownames(var = "gene_symbol") %>%
  t() %>% scale() %>% t() %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "gene_symbol") %>%
  transpose_df(genelist = test_genelist)

Maynard_subset <- Maynard_logCPM %>% 
  column_to_rownames(var = "gene_symbol") %>%
  t() %>% scale() %>% t() %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "gene_symbol") %>%
  transpose_df(genelist = test_genelist)

He_cor <- cor(He_subset, He_subset, method = "pearson")
He_cor[is.na(He_cor)] <- 0
Maynard_cor <- cor(Maynard_subset, Maynard_subset, method = "pearson")

mantel(He_cor, Maynard_cor, method = "spearman", permutations = 1)
# Mantel statistic r: 0.5882