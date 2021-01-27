#Script to scale and export MTG scRNA-seq data from Allen

library(tidyverse)
library(tidyr)
library(data.table)
library(moments)
library(here)
library(magrittr)
library(edgeR)

# Functions
# Read in data and metadata

scMetadata <- fread(here("Data", "Allen", "singlecellMetadata.csv"), header = T) %>%
  dplyr::select(sample_name, class_label, subclass_label, region_label, cortical_layer_label, outlier_call) %>%
  as_tibble()
scMatrix <- fread(here("Data", "Allen", "singlecellMatrix.csv"), header = T)

normalize_AIBS <- function(AIBS_countMatrix, AIBS_metadata) {
  metadata <- AIBS_metadata
  matrix <- AIBS_countMatrix
  
  matrix_logCPM <- matrix %>%
    column_to_rownames(var = "sample_name") %>%
    t() %>%
    cpm(log = T) %>%
    t() %>% as.data.frame() %>%
    rownames_to_column(var = "sample_name")
  
  matrix_logCPM$class_label <- metadata$class_label
  matrix_logCPM$region_label <- metadata$region_label
  matrix_logCPM$outlier_call <- metadata$outlier_call
  matrix_logCPM$cortical_layer_label <- metadata$cortical_layer_label
  
  matrix_logCPM %<>% filter(outlier_call == FALSE)
  return(matrix_logCPM)
}

select_AIBS_region <- function(normalized_AIBS_df, region) {
  
  matrix <- normalized_AIBS_df
  matrix %<>% filter(region_label == region) %>%
    gather("gene", "expression_value", -c("sample_name", "class_label", "region_label", "cortical_layer_label"))
  
  return(matrix)
}

cellType_group <- function(matrix_region, cellType) {
  cellTypeSpecific_matrix <- matrix_region %>%
    filter(class_label == cellType) %>%
    select(-class_label, -region_label, -sample_name) %>%
    group_by(gene, cortical_layer_label) %>%
    summarise(mean_expression = mean(expression_value)) %>%
    rename(gene_symbol = gene) %>%
    add_column(class_label = cellType) %>%
    select(gene_symbol, class_label, cortical_layer_label, mean_expression)
  return(cellTypeSpecific_matrix)
}

# Run to obtain 

scMatrix_logCPM <- normalize_AIBS(scMatrix, scMetadata)
scMatrix_logCPM_MTG <- select_AIBS_region(scMatrix_logCPM, "MTG")


MTG_matrix_scaled <- list()
for (i in c("GABAergic", "Glutamatergic", "Non-neuronal")) {
  MTG_matrix_scaled[[i]] <- cellType_group(scMatrix_logCPM_MTG, i)
}

MTG_matrix_scaled <- rbind(MTG_matrix_scaled$GABAergic, MTG_matrix_scaled$Glutamatergic, MTG_matrix_scaled$`Non-neuronal`)

save(MTG_df, file = here("Data", "Allen", "MTG_df_01_21.Rdata"))
write.csv(MTG_df, file = here("Data", "Allen", "MTG_df_01_21.Rdata"))
