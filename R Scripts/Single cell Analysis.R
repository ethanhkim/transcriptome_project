library(tidyverse)
library(here)
library(data.table)
library(magrittr)
library(dplyr)
library(tidyr)


# Download Allen sc-RNAseq matrix and metadata for CAMH SCC

singlecellMetadata_dest <- here("Data", "Allen", "singlecellMetadata.csv")
singlecellMetadata_url <- "https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/metadata.csv"
allen_singlecellMetadata <- download.file(singlecellMetadata_url, singlecellMetadata_dest)

singlecellMatrix_dest <- here("Data", "Allen", "singlecellMatrix.csv")
singlecellMatrix_url <- "https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/matrix.csv"
allen_singlecellMatrix <- download.file(singlecellMatrix_url, singlecellMatrix_dest)

<<<<<<< HEAD
=======
<<<<<<< HEAD

=======
>>>>>>> 1936333483675149dafcbe0a15a472735421bf81
>>>>>>> 2a65a83329b2bf7d1a56ddb0a32ab14e1c7bded5
# Function to read in Allen sc-RNAseq matrix and associated metadata separate the data into regions
separate_by_region <- function(region) {
  
  metadata <- fread(here("Data", "Allen", "singlecellMetadata.csv"), header = T) %>%
    dplyr::select(sample_name, class_label, subclass_label, region_label, cortical_layer_label, outlier_call) %>%
    as_tibble()
  matrix <- fread(here("Data", "Allen", "singlecellMatrix.csv"), header = T) %>%
    as_tibble()
  
  matrix$class_label <- metadata$class_label
  matrix$region_label <- metadata$region_label
  matrix$cortical_layer_label <- metadata$cortical_layer_label
  matrix$outlier_call <- metadata$outlier_call
  
  matrix %<>%
    dplyr::select(sample_name, class_label, region_label, cortical_layer_label, everything()) %>%
    filter(outlier_call == FALSE) %>%
    filter(region_label == region)
  metadata <- matrix %>%
    dplyr::select(sample_name:cortical_layer_label)
  matrix %<>%
    dplyr::select(-class_label, -region_label, -cortical_layer_label) %>%
    column_to_rownames(var = "sample_name") %>%
    scale() %>%
    as.data.frame() %>%
    rownames_to_column(var = "sample_name") %>%
    merge(metadata, by = "sample_name") %>%
    gather("gene", "expression_value", -c("sample_name", "class_label", "region_label", "cortical_layer_label")) %>%
    dplyr::select(sample_name, gene, class_label, region_label, cortical_layer_label, expression_value) %>%
    group_by(gene, class_label, region_label, cortical_layer_label) %>%
    summarize(median_exp_value = median(expression_value))
  
}

<<<<<<< HEAD
=======
<<<<<<< HEAD
#Separate by region
=======
<<<<<<< HEAD
# Separate into region-specific data in long-form
=======
>>>>>>> 0da6db4d7849688623b51264de08d0f96a5c98c6

>>>>>>> 1936333483675149dafcbe0a15a472735421bf81
>>>>>>> 2a65a83329b2bf7d1a56ddb0a32ab14e1c7bded5
A1C_matrix <- separate_by_region("A1C")
MTG_matrix <- separate_by_region("MTG")
V1C_matrix <- separate_by_region("V1C")
CgG_matrix <- separate_by_region("CgG")
M1lm_matrix <- separate_by_region("M1lm")
M1ul_matrix <- separate_by_region("M1ul")
S1ul_matrix <- separate_by_region("S1ul")
S1lm_matrix <- separate_by_region("S1lm")

<<<<<<< HEAD

#Write matrices to csv

# Save region-specific long-form data in CSV for webapp

=======
<<<<<<< HEAD
#Write matrices to csv
=======
<<<<<<< HEAD
# Save region-specific long-form data in CSV for webapp
=======

>>>>>>> 0da6db4d7849688623b51264de08d0f96a5c98c6

>>>>>>> 1936333483675149dafcbe0a15a472735421bf81
>>>>>>> 2a65a83329b2bf7d1a56ddb0a32ab14e1c7bded5
write.csv(A1C_matrix, here("Data", "Allen", "A1C_matrix.csv"))
write.csv(MTG_matrix, here("Data", "Allen", "MTG_matrix.csv"))
write.csv(V1C_matrix, here("Data", "Allen", "V1C_matrix.csv"))
write.csv(CgG_matrix, here("Data", "Allen", "CgG_matrix.csv"))
write.csv(M1lm_matrix, here("Data", "Allen", "M1lm_matrix.csv"))
write.csv(M1ul_matrix, here("Data", "Allen", "M1ul_matrix.csv"))
write.csv(S1ul_matrix, here("Data", "Allen", "S1ul_matrix.csv"))
write.csv(S1lm_matrix, here("Data", "Allen", "S1lm_matrix.csv"))
<<<<<<< HEAD
=======


<<<<<<< HEAD
=======


>>>>>>> 0da6db4d7849688623b51264de08d0f96a5c98c6

>>>>>>> 2a65a83329b2bf7d1a56ddb0a32ab14e1c7bded5
