#Script to scale and export MTG scRNA-seq data from Allen

library(tidyverse)
library(tidyr)
library(data.table)
library(moments)
library(here)

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

region_specific_scRNA_df <- function(region) {
  
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
  
  allen_MTG_matrix <- matrix %>%
    filter(region_label == region) %>%
    gather("gene", "expression_value", -c("sample_name", "class_label", "region_label", "cortical_layer_label"))
}



GABA_MTG <- allen_MTG_matrix %>% 
  filter(class_label == "GABAergic") %>%
  select(-class_label, -region_label, -sample_name) %>%
  group_by(gene, cortical_layer_label) %>%
  summarise(median_expression = median(expression_value)) %>%
  spread(cortical_layer_label, median_expression) %>%
  rename(gene_symbol = gene) %>%
  ungroup() %>%
  column_to_rownames(var = "gene_symbol") %>% 
  t() %>%
  scale() %>%
  t() %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "gene_symbol") %>%
  add_column(class_label = "GABAergic") %>%
  select(gene_symbol, class_label, everything())


GLUT_MTG <- allen_MTG_matrix %>% 
  filter(class_label == "Glutamatergic") %>%
  select(-class_label, -region_label, -sample_name) %>%
  group_by(gene, cortical_layer_label) %>%
  summarise(median_expression = median(expression_value)) %>%
  spread(cortical_layer_label, median_expression) %>%
  rename(gene_symbol = gene) %>%
  ungroup() %>%
  column_to_rownames(var = "gene_symbol") %>% 
  t() %>%
  scale() %>%
  t() %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "gene_symbol") %>%
  add_column(class_label = "Glutamatergic") %>%
  select(gene_symbol, class_label, everything())

NON_MTG <- allen_MTG_matrix %>% 
  filter(class_label == "Non-neuronal") %>%
  select(-class_label, -region_label, -sample_name) %>%
  group_by(gene, cortical_layer_label) %>%
  summarise(median_expression = median(expression_value)) %>%
  spread(cortical_layer_label, median_expression) %>%
  rename(gene_symbol = gene) %>%
  ungroup() %>%
  column_to_rownames(var = "gene_symbol") %>% 
  t() %>%
  scale() %>%
  t() %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "gene_symbol") %>%
  add_column(class_label = "Non-neuronal") %>%
  select(gene_symbol, class_label, everything())

MTG_matrix_scaled <- rbind(GABA_MTG, GLUT_MTG, NON_MTG)

save(MTG_matrix_scaled, file = here("Data", "Allen", "MTG_matrix_scaled.Rdata"))



