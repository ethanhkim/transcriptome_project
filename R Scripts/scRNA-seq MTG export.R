#Script to scale and export MTG scRNA-seq data from Allen

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}


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
  filter(outlier_call == FALSE) %>%
  filter(region_label == "MTG")

metadata <- matrix %>%
  dplyr::select(sample_name:cortical_layer_label)

allen_matrix_MTG <- matrix %>%
  gather("gene", "expression_value", -c("sample_name", "class_label", "region_label", "cortical_layer_label"))

set.seed(1)

test <- allen_matrix_MTG %>%
  group_by(class_label) %>%
  group_by(gene) %>%
  mutate(expression_value = scale_this(expression_value)) %>%
  group_by(gene, class_label, region_label, cortical_layer_label) %>%
  summarize(median_exp_value = median(expression_value, na.rm = FALSE))

test_genelist <- c("RELN", "NDNF", "CHRNA7", "CNR1", "CXCL14", "DISC1")

test_data <- test %>%
  filter(gene %in% test_genelist)

test %>%
  filter(class_label == "Non-neuronal") %>%
  filter(gene %in% test_genelist) %>%
  ggplot(mapping = aes(x = cortical_layer_label, y = gene, fill = median_exp_value)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu", limits = c(-1,1)*max(abs(test$median_exp_value))) +
  scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand=c(0,0))

group_by(class_label) %>%
  group_by(gene) %>%
  scale_this()
dplyr::select(-class_label, -region_label, -cortical_layer_label) %>%
  column_to_rownames(var = "sample_name") %>%
  group_by(class_label) %>%
  group_by(gene) %>%
  scale_this()
scale() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample_name") %>%
  merge(metadata, by = "sample_name") %>%
  gather("gene", "expression_value", -c("sample_name", "class_label", "region_label", "cortical_layer_label")) %>%
  dplyr::select(sample_name, gene, class_label, region_label, cortical_layer_label, expression_value) %>%
  group_by(gene, class_label, region_label, cortical_layer_label) %>%
  summarize(median_exp_value = median(expression_value, na.rm = FALSE))



