#Script to scale and export MTG scRNA-seq data from Allen

library(tidyverse)
library(tidyr)
library(data.table)
library(moments)
library(here)
library(magrittr)

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

allen_MTG_matrix <- region_specific_scRNA_df("MTG")

scale_scRNA_region <- function(scRNA_data, cellType) {
  scaled_data <- scRNA_data %>%
    filter(class_label == cellType) %>%
    select(-class_label, -region_label, -sample_name) %>%
    group_by(gene) %>%
    mutate(expression_log = log(expression_value + 1)) %>%
    ungroup() %>%
    group_by(gene, cortical_layer_label) %>%
    summarise(median_expression = median(expression_log)) %>%
    mutate(median_expression_log_scaled = scale_this(median_expression)) %>%
    rename(gene_symbol = gene) %>%
    add_column(class_label = cellType) %>%
    select(gene_symbol, class_label, cortical_layer_label, median_expression)
}

MTG_matrix_scaled <- list()
for (i in c("GABAergic", "Glutamatergic", "Non-neuronal")) {
  MTG_matrix_scaled[[i]] <- scale_scRNA_region(allen_MTG_matrix, i)
}

test <- MTG_matrix_scaled
MTG_matrix_scaled %<>% flatten(.id = "cell_type")

GABA_MTG <- scale_scRNA_region(allen_MTG_matrix, "GABAergic")
GLUT_MTG <- scale_scRNA_region(allen_MTG_matrix, "Glutamatergic")
GABA_MTG <- scale_scRNA_region(allen_MTG_matrix, "GABAergic")

MTG_matrix_scaled <- rbind(MTG_matrix_scaled$GABAergic, MTG_matrix_scaled$Glutamatergic, MTG_matrix_scaled$`Non-neuronal`)

save(MTG_matrix_scaled, file = here("Data", "Allen", "MTG_matrix_scaled.Rdata"))



