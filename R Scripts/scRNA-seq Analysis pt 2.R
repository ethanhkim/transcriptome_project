library(moments)

metadata <- fread(here("Data", "Allen", "singlecellMetadata.csv"), header = T) %>%
  dplyr::select(sample_name, class_label, subclass_label, region_label, cortical_layer_label, outlier_call) %>%
  as_tibble()
matrix <- fread(here("Data", "Allen", "singlecellMatrix.csv"), header = T, integer64 = "numeric") %>%
  as_tibble()


matrix$class_label <- metadata$class_label
matrix$region_label <- metadata$region_label
matrix$cortical_layer_label <- metadata$cortical_layer_label
matrix$outlier_call <- metadata$outlier_call

matrix %<>%
  dplyr::select(sample_name, class_label, region_label, cortical_layer_label, everything()) %>%
  filter(outlier_call == FALSE)
metadata <- matrix %>%
  dplyr::select(sample_name:cortical_layer_label)

allen_matrix_long <- matrix %>%
  dplyr::select(-class_label, -region_label, -cortical_layer_label) %>%
  column_to_rownames(var = "sample_name") %>%
  scale() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample_name") %>%
  merge(metadata, by = "sample_name") %>%
  gather("gene", "expression_value", -c("sample_name", "class_label", "region_label", "cortical_layer_label")) %>%
  dplyr::select(sample_name, gene, class_label, region_label, cortical_layer_label, expression_value) %>%
  group_by(gene, class_label, region_label, cortical_layer_label) %>%
  summarize(median_exp_value = median(expression_value, na.rm = FALSE))

testmatrix[is.na(testmatrix)] <- 0
skewness(testmatrix$median_exp_value)
# severe skew - coefficient of 124.2559

testmatrix$median_exp_value <- log10(testmatrix$median_exp_value)


