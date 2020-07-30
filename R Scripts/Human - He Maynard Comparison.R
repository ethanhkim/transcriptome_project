# More comparions

library(tidyverse)
library(magrittr)

He_Maynard_common_gene_list <- intersect(He_DS1_averaged_by_layer$gene_symbol, Maynard_dataset_average$gene_symbol)
testlist <- c("GAD1", "CCK", "GRIN1")

He_cormatrix <- He_DS1_averaged_by_layer %>%
  filter(gene_symbol %in% testlist) %>%
  arrange(gene_symbol) %>%
  column_to_rownames(var = "gene_symbol") %>%
  t()

Maynard_cormatrix <- Maynard_dataset_average %>%
  select(gene_symbol:WM) %>%
  filter(gene_symbol %in% testlist) %>%
  column_to_rownames(var = "gene_symbol") %>%
  t()

He_Maynard_cormatrix <- cor(He_cormatrix, Maynard_cormatrix, method = "pearson") 
He_Maynard_diag <- diag(He_Maynard_cormatrix, names = TRUE)

He_Maynard_vector <- format(round(He_Maynard_diag, 2), nsmall = 2) %>%
  as_tibble(.name_repair = "unique") %>%
  pull(var = 1) %>%
  as.numeric() %>%
  mean() 

hist(He_Maynard_diag)

test <- ecdf(He_Maynard_diag)
test(0.79)

quantile(He_Maynard_diag, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)

He_Maynard_vector <- He_Maynard_cormatrix %>%
  pull(var = 1) %>%
  as.numeric()
He_Maynard_vector <- format(round(He_Maynard_vector, 2), nsmall = 2) %>%
  as_tibble()

