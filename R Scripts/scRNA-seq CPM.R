## Allen SMART-seq Multiple Regions CPM normalization ##

# Load required libraries
library(tidyverse)
library(tidyr)
library(data.table)
library(moments)
library(here)
library(magrittr)
library(edgeR)
library(conflicted)

# Set conflicts
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("transpose", "data.table")

# Download Allen sc-RNAseq matrix and metadata for CAMH SCC
singlecellMetadata_dest <- here("Data", "raw_data", "Allen", "singlecellMetadata.csv")
singlecellMetadata_url <- "https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/metadata.csv"
allen_singlecellMetadata <- download.file(singlecellMetadata_url, singlecellMetadata_dest)

singlecellMatrix_dest <- here("Data", "raw_data", "Allen", "singlecellMatrix.csv")
singlecellMatrix_url <- "https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/matrix.csv"
allen_singlecellMatrix <- download.file(singlecellMatrix_url, singlecellMatrix_dest)

rm(singlecellMetadata_dest, singlecellMetadata_url,
   singlecellMatrix_dest, singlecellMatrix_url,
   allen_singlecellMatrix, allen_singlecellMetadata)

# Functions to use #

# Read in data and metadata
AIBS_metadata <- fread(here("Data", "raw_data", "Allen", 
                            "singlecellMetadata.csv"), header = T) %>%
  select(sample_name, class_label, subclass_label, region_label, 
         cortical_layer_label, outlier_call) %>%
  as_tibble()
AIBS_matrix <- fread(here("Data", "raw_data", "Allen",
                          "singlecellMatrix.csv"), header = T)

normalize_AIBS <- function(AIBS_countMatrix, AIBS_metadata) {
  metadata <- AIBS_metadata
  matrix <- AIBS_countMatrix
  
  matrix$sample_name <- NULL
  # Get gene symbols from columns
  matrix_colnames <- colnames(matrix)
  
  matrix_CPM <- matrix %>%
    transpose() %>%
    cpm() %>% as.data.frame() %>%
    transpose() %>% as.data.frame()
  
  colnames(matrix_CPM) <- matrix_colnames
  
  matrix_CPM$class_label <- metadata$class_label
  matrix_CPM$region_label <- metadata$region_label
  matrix_CPM$outlier_call <- metadata$outlier_call
  matrix_CPM$cortical_layer_label <- metadata$cortical_layer_label
  matrix_CPM$sample_name <- metadata$sample_name
  
  matrix_CPM %<>% filter(outlier_call == FALSE)
  return(matrix_CPM)
}

select_AIBS_region <- function(AIBS_CPM_df, region) {
  
  matrix <- AIBS_CPM_df
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
    mutate(log_mean_expression = log2(mean_expression + 1)) %>%
    rename(gene_symbol = gene) %>%
    add_column(class_label = cellType) %>%
    select(gene_symbol, class_label, cortical_layer_label, log_mean_expression)
  return(cellTypeSpecific_matrix)
}

# Run to obtain 
AIBS_CPM <- normalize_AIBS(AIBS_matrix, AIBS_metadata)
# Create region-specific dataframes
create_region_df <- function(region_df) {
  
  df_list <- list()
  for (i in c("GABAergic", "Glutamatergic", "Non-neuronal")) {
    df_list[[i]] <- cellType_group(region_df, i)
  }
  
  df <- rbind(df_list$GABAergic, df_list$Glutamatergic,
              df_list$`Non-neuronal`)
  
  return(df)
}

# Middle temporal gyrus
MTG <- select_AIBS_region(AIBS_CPM, "MTG")
MTG_df <- create_region_df(MTG) 
MTG_df$region_label <- "MTG"
rm(MTG)

# Primary motor cortex
M1lm <- select_AIBS_region(AIBS_CPM, "M1lm")
M1lm_df <- create_region_df(M1lm)
M1lm_df$region_label <- "M1lm"
M1lm_df %<>% 
rm(M1lm)
M1ul <- select_AIBS_region(AIBS_CPM, "M1ul")
M1ul_df <- create_region_df(M1ul)
M1ul_df$region_label <- "M1ul"
rm(M1ul)

# Primary sensory cortex
S1lm <- select_AIBS_region(AIBS_CPM, "S1lm")
S1lm_df <- create_region_df(S1lm)
S1lm_df$region_label <- "S1lm"
rm(S1lm)
S1ul <- select_AIBS_region(AIBS_CPM, "S1ul")
S1ul_df <- create_region_df(S1ul)
S1ul_df$region_label <- "S1ul"
rm(S1ul)

# Primary auditory cortex
A1C <- select_AIBS_region(AIBS_CPM, "A1C")
A1C_df <- create_region_df(A1C)
A1C_df$region_label <- "A1C"
rm(A1C)

# Primary visual cortex
V1C <- select_AIBS_region(AIBS_CPM, "V1C")
V1C_df <- create_region_df(V1C)
V1C_df$region_label <- "V1C"
rm(V1C)

# Cingulate gyrus
CgG <- select_AIBS_region(AIBS_CPM, "CgG")
CgG_df <- create_region_df(CgG)
CgG_df$region_label <- "CgG"
CgG_df %<>%
  group_by(gene_symbol, cortical_layer_label) %>%
  summarize(mean_exp = mean(log_mean_expression)) %>%
  spread(cortical_layer_label, mean_exp) %>%
  mutate(L5 = (L5a+L5b)/2) %>%
  add_column(L3 = 0, L4 = 0) %>%
  select(gene_symbol, L1, L2, L3, L4, L5, L6)
rm(CgG)

# Remove CPM dataframe
rm(AIBS_CPM)
# Create logCPM dataframe
AIBS_logCPM <- rbind(A1C_df, CgG_df, M1lm_df, M1ul_df, 
                     MTG_df, S1lm_df, S1ul_df, V1C_df)

# Save logCPM data
save(AIBS_logCPM, file = "/external/mgmt3/genome/scratch/Neuroinformatics/ekim/transcriptome_project/Data/processed_data/AIBS_logCPM.RData")
