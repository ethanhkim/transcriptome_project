## Allen SMART-seq Multiple Regions CPM normalization ##

# Load required libraries
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
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

# Check if working at SCC or Local
if (getwd() == "/external/mgmt3/genome/scratch/Neuroinformatics/ekim/transcriptome_project") {
  # Check if Allen count matrix and metadata exist
  if (!file.exists(here("Data", "raw_data", "Allen",
                        "singlecellMatrix.csv")) == TRUE) {
    singlecellMatrix_dest <- here("Data", "raw_data", "Allen", "singlecellMatrix.csv")
    singlecellMatrix_url <- "https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/matrix.csv"
    allen_singlecellMatrix <- download.file(singlecellMatrix_url, singlecellMatrix_dest)
    
    rm(singlecellMatrix_dest, singlecellMatrix_url, allen_singlecellMatrix)
  }
  
  if (!file.exists(here("Data", "raw_data", "Allen",
                        "singlecellMetadata.csv")) == TRUE) {
    # Download Allen sc-RNAseq matrix and metadata for CAMH SCC
    singlecellMetadata_dest <- here("Data", "raw_data", "Allen", "singlecellMetadata.csv")
    singlecellMetadata_url <- "https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/metadata.csv"
    allen_singlecellMetadata <- download.file(singlecellMetadata_url, singlecellMetadata_dest)
    
    rm(singlecellMetadata_dest, singlecellMetadata_url, allen_singlecellMetadata)
  }
}


# Functions to use #

# Create AIBS data with merged metadata and remove
# outliers
clean_AIBS_data <- function(AIBS_countMatrix, AIBS_metadata) {
  metadata <- AIBS_metadata
  matrix <- AIBS_countMatrix
  
  matrix$outlier_call <- metadata$outlier_call
  matrix$class_label <- metadata$class_label
  matrix$region_label <- metadata$region_label
  matrix$cortical_layer_label <- metadata$cortical_layer_label
  matrix$external_donor_name_label <- metadata$external_donor_name_label
  
  matrix %<>% filter(outlier_call == FALSE) %>%
    select(-outlier_call) 
  return(matrix)
}

# Summarize gene counts by layer and cell-type
sum_gene_count <- function(region_separated_df) {
  
  df <- region_separated_df
  
  df %<>%
    select(sample_name, class_label, cortical_layer_label,
           external_donor_name_label, everything()) %>%
    # Lengthen dataframe to:
    # gene_symbol | count | everything()
    data.table::melt(id.vars = c("sample_name", "class_label", 
                                 "cortical_layer_label",
                                 "external_donor_name_label"),
                     variable.factor = FALSE,
                     variable.name = "gene_symbol",
                     value.name = "count") %>%
    # Group by layer, cell type (class_label) and gene 
    group_by(cortical_layer_label, class_label,
             gene_symbol) %>%
    # Sum up all gene counts at a layer and cell type level
    summarize(countSum = sum(count)) %>%
    # Widen dataframe to gene (row) by layer (column)
    spread(cortical_layer_label, countSum)
  
  return(df)
}

# CPM normalize summarized count matrix
cpm_normalize <- function(summed_count_matrix) {
  
  df <- summed_count_matrix
  # Save gene symbol vector
  df_gene_symbol <- df$gene_symbol
  # Save class_label vector
  df_class_label <- df$class_label
  # Remove to create numerical matrix for cpm()
  df$gene_symbol <- NULL
  df$class_label <- NULL
  # logCPM normalize and add back in gene and cell type labels
  df %<>% cpm(log = T) %>% as.data.frame() %>%
    add_column(gene_symbol = df_gene_symbol,
               class_label = df_class_label) %>%
    select(class_label, gene_symbol, everything())
  
  return(df)
}

# Read in data and metadata
AIBS_metadata <- fread(here("Data", "raw_data", "Allen", 
                            "singlecellMetadata.csv"), header = T) %>%
  select(sample_name, external_donor_name_label, class_label, subclass_label, region_label, 
         cortical_layer_label, outlier_call) %>%
  as_tibble()
AIBS_matrix <- fread(here("Data", "raw_data", "Allen",
                          "singlecellMatrix.csv"), header = T)

# Create cleaned dataframe and remove source data
AIBS_cleaned_df <- clean_AIBS_data(AIBS_matrix, AIBS_metadata)
rm(AIBS_matrix, AIBS_metadata)

# Separate by region (MTG)
MTG <- AIBS_cleaned_df %>%
  filter(region_label == "MTG") %>%
  select(-region_label)
rm(AIBS_cleaned_df)

MTG_sum_count <- sum_gene_count(MTG)
MTG_logCPM_dataset <- cpm_normalize(MTG_sum_count) 
MTG_logCPM_dataset %<>% add_column(WM = NA)

# Save logCPM data
save(MTG_logCPM_dataset, file = here("Data", "processed_data", "MTG_logCPM_dataset.Rdata"))
write.csv(MTG_logCPM_dataset, file = here("Data", "processed_data", "MTG_logCPM_dataset.csv"))
