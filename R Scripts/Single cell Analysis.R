library(tidyverse)
library(here)
library(data.table)
library(magrittr)
library(dplyr)
library(tidyr)
library(here)

allen_singlecellMetadata <- fread(here("Data", "Allen", "metadata.csv"), header = T) %>%
  select(sample_name, class_label, subclass_label, region_label, cortical_layer_label, outlier_call) %>%
  as_tibble()
allen_singlecellMatrix <- fread(here("Data", "Allen", "matrix.csv"), header = T) %>%
  as_tibble()

allen_singlecellMatrix$class_label <- allen_singlecellMetadata$class_label
allen_singlecellMatrix$region_label <- allen_singlecellMetadata$region_label
allen_singlecellMatrix$cortical_layer_label <- allen_singlecellMetadata$cortical_layer_label
allen_singlecellMatrix$outlier_call <- allen_singlecellMetadata$outlier_call

allen_singlecellMatrix %<>%
  select(sample_name, class_label, region_label, cortical_layer_label, everything()) %>%
  filter(outlier_call == FALSE) 

separate_by_region <- function(allenMatrix, region) {
  
  metadata <- fread(here("Data", "Allen", "metadata.csv"), header = T) %>%
    select(sample_name, class_label, subclass_label, region_label, cortical_layer_label, outlier_call) %>%
    as_tibble()
  matrix <- fread(here("Data", "Allen", "matrix.csv"), header = T) %>%
    as_tibble()
  
  matrix$class_label <- metadata$class_label
  matrix$region_label <- metadata$region_label
  matrix$cortical_layer_label <- metadata$cortical_layer_label
  matrix$outlier_call <- metadata$outlier_call
  
  matrix %<>%
    select(sample_name, class_label, region_label, cortical_layer_label, everything()) %>%
    filter(outlier_call == FALSE) 
  
  matrix %<>%
    filter(region_label == region, .keep_all = TRUE)
  metadata <- matrix %>%
    select(sample_name:cortical_layer_label)
  matrix %<>%
    select(-class_label, -subclass_label, -region_label, -cortical_layer_label) %>%
    column_to_rownames(var = "sample_name") %>%
    scale() %>%
    as.data.frame() %>%
    rownames_to_column(var = "sample_name") %>%
    merge(metadata, by = "sample_name") %>%
    gather("gene", "expression_value", -c("sample_name", "class_label", "region_label", "cortical_layer_label")) %>%
    select(sample_name, class_label, subclass_label, region_label, cortical_layer_label, everything())

}

testmatrix <- head(allen_singlecellMatrix, n = 20)

allen_singlecellMatrix_A1C <- separate_by_region(allen_singlecellMatrix, "A1C")

allen_singlecellMatrix_A1C <- allen_singlecellMatrix %>%
  filter(region_label == "A1C", .keep_all = TRUE)
A1C_metadata <- allen_singlecellMatrix_A1C %>%
  select(sample_name:cortical_layer_label)
allen_singlecellMatrix_A1C %<>%
  select(-class_label, )


allen_singlecellMatrix <- rbind(allen_singlecellMatrix_A1C, allen_singlecellMatrix_CgG, allen_singlecellMatrix_M1lm, allen_singlecellMatrix_M1ul, 
                                allen_singlecellMatrix_MTG, allen_singlecellMatrix_S1lm, allen_singlecellMatrix_S1ul, allen_singlecellMatrix_V1C)
  
