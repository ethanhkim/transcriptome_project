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


# Function to read in Allen sc-RNAseq matrix and associated metadata separate the data into regions
separate_by_region <- function(region) {
  
  metadata <- fread("/Users/ethankim/Desktop/transcriptome_project/Data/Allen/metadata.csv", header = T) %>%
    dplyr::select(sample_name, class_label, subclass_label, region_label, cortical_layer_label, outlier_call) %>%
    as_tibble()
  matrix <- fread("/Users/ethankim/Desktop/transcriptome_project/Data/Allen/matrix.csv", header = T) %>%
    as_tibble()
  
  matrix$class_label <- metadata$class_label
  matrix$region_label <- metadata$region_label
  matrix$cortical_layer_label <- metadata$cortical_layer_label
  matrix$outlier_call <- metadata$outlier_call
  
  matrix %<>%
    dplyr::select(sample_name, class_label, region_label, cortical_layer_label, everything()) %>%
    filter(outlier_call == FALSE) %>%
    filter(region_label == region)
  
  return(matrix)
  
}

metadata <- fread("/Users/ethankim/Desktop/transcriptome_project/Data/Allen/metadata.csv", header = T) %>%
  dplyr::select(sample_name, class_label, subclass_label, region_label, cortical_layer_label, outlier_call) %>%
  as_tibble()
MTG_matrix <- fread("/Users/ethankim/Desktop/transcriptome_project/Data/Allen/matrix.csv", header = T) %>%
  as_tibble()


#Separate by region
MTG_matrix <- separate_by_region("MTG")





#read in CSV
M1ul_matrix <- read.csv(here("Data", "Allen", "M1ul_matrix.csv")) %>%
  select(-X) %>%
  unite("gene_id", gene:cortical_layer_label, sep = ",") %>%
  column_to_rownames(var = "gene_id")
  

M1lm_matrix <- read.csv(here("Data", "Allen", "M1lm_matrix.csv")) %>%
  select(-X) %>%
  unite("gene_id", gene:cortical_layer_label, sep = ",") %>%
  column_to_rownames(var = "gene_id")

M1_matrix <- M1ul_matrix %>%
  add_column(median_exp_val_M1lm = M1lm_matrix$median_exp_value)
M1_matrix$mean_exp <- rowMeans(M1_matrix[1:2], na.rm = T) 
M1_matrix %<>%
  select(mean_exp) %>%
  rownames_to_column("labels") %>%
  separate("labels", into = c("gene", "class_label", "region_label", "cortical_layer_label"), sep = ",") %>%
  mutate(region_label = gsub("1ul", "1", region_label)) %>%
  rename("expression" = mean_exp)


#read in CSV
S1ul_matrix <- read.csv(here("Data", "Allen", "S1ul_matrix.csv")) %>%
  select(-X) %>%
  unite("gene_id", gene:cortical_layer_label, sep = ",") %>%
  column_to_rownames(var = "gene_id")


S1lm_matrix <- read.csv(here("Data", "Allen", "S1lm_matrix.csv")) %>%
  select(-X) %>%
  unite("gene_id", gene:cortical_layer_label, sep = ",") %>%
  column_to_rownames(var = "gene_id")

S1_matrix <- S1ul_matrix %>%
  add_column(median_exp_val_S1lm = S1lm_matrix$median_exp_value)
S1_matrix$mean_exp <- rowMeans(S1_matrix[1:2], na.rm = T) 
S1_matrix %<>%
  select(mean_exp) %>%
  rownames_to_column("labels") %>%
  separate("labels", into = c("gene", "class_label", "region_label", "cortical_layer_label"), sep = ",") %>%
  mutate(region_label = gsub("1ul", "1", region_label)) %>%
  mutate_at(vars(mean_exp), ~replace(., is.nan(.), NA)) %>%
  rename("expression" = mean_exp)

write.csv(M1_matrix, here("Data", "Allen", "M1_matrix.csv")) 
write.csv(S1_matrix, here("Data", "Allen", "S1_matrix.csv"))
  

A1C_matrix <- read.csv(here("Data", "Allen", "A1C_matrix.csv")) %>%
  select(-X) %>%
  rename("expression" = median_exp_value)
MTG_matrix <- read.csv(here("Data", "Allen", "MTG_matrix.csv")) %>%
  select(-X) %>%
  rename("expression" = median_exp_value)
V1C_matrix <- read.csv(here("Data", "Allen", "V1C_matrix.csv")) %>%
  select(-X) %>%
  rename("expression" = median_exp_value)
CgG_matrix <- read.csv(here("Data", "Allen", "CgG_matrix.csv")) %>%
  select(-X) %>%
  rename("expression" = median_exp_value)


merged_allen_scRNA_data_long <- rbind(A1C_matrix, MTG_matrix, V1C_matrix, CgG_matrix, S1_matrix, M1_matrix)
merged_allen_scRNA_data_matrix <- merged_allen_scRNA_data_long %>%
  pivot_wider(
  names_from = cortical_layer_label,
  values_from = expression
) %>%
  select(gene, class_label, region_label, L1, L2, L3, L4, L4ab, L5, L5a, L5b, L6, L6a, L6b, WM) %>%
  rename(gene_symbol = gene) %>%
  rowwise() %>%
  mutate(L4 = mean(c(L4, L4ab), na.rm = T)) %>%
  mutate(L5 = mean(c(L5, L5a, L5b), na.rm = T)) %>%
  mutate(L6 = mean(c(L6, L6a, L6b), na.rm = T)) %>%
  select(gene_symbol, class_label, region_label, L1, L2, L3, L4, L5, L6, WM) %>%
  mutate_at(vars(contains("L")), ~replace(., is.nan(.), NA)) %>%
  rename(Layer_1 = L1, Layer_2 = L2, Layer_3 = L3, Layer_4 = L4, Layer_5 = L5, Layer_6 = L6)

write.csv(merged_allen_scRNA_data_long, here("Data", "Allen", "merged_Allen_scRNA_data_long.csv"))
write.csv(merged_allen_scRNA_data_matrix, here("Data", "Allen", "merged_Allen_scRNA_data_matrix.csv"))
save(merged_allen_scRNA_data_matrix, file = here("Data", "Allen", "merged_Allen_scRNA_data_matrix.Rdata"))
save(MTG_matrix, file = here("Data", "Allen", "MTG_matrix.Rdata"))





