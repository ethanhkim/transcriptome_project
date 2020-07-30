He_DS1_Human_values <- He_DS1_Human %>%
  column_to_rownames(var = "gene_symbol") %>%
  select(-"row_index")

He_DS2_Human_values <- He_DS2_Human %>%
  column_to_rownames(var = "gene_symbol") %>%
  lapply(function(x) as.numeric(as.character(x))) %>%
  as.data.frame()
rownames(He_DS1_Human_values) <- He_DS2_Human$gene_symbol

#Transpose table such that table represents the 17 separate cuts for each gene
He_DS1_Human_values_transposed <- He_DS1_Human_values %>%
  t() %>%
  as_tibble() %>%
  replace(is.na(.), 0)

He_DS2_Human_values_transposed <- He_DS2_Human_values %>%
  t() %>%
  as_tibble() %>%
  replace(is.na(.), 0)
  

# Create correlation matrices
He_DS1_cormatrix <- stats::cor(He_DS1_Human_values_transposed, He_DS1_Human_values_transposed, method = "pearson") %>%
  replace(is.na(.), 0)
He_DS2_cormatrix <- stats::cor(He_DS2_Human_values_transposed, He_DS2_Human_values_transposed, method = "pearson") %>%
  replace(is.na(.), 0)


He_manteltest <- mantel(He_DS1_cormatrix, He_DS2_cormatrix, method = "spearman", permutations = 1) %>%
  print()

Mantel_test_result_statistic <- results$statistic

