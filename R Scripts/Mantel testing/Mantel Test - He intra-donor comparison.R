## He et al Normalization ##

## Load in libraries ----
library(edgeR)
library(data.table)
library(dplyr)
library(purrr)
library(stringr)
library(here)

## Functions to use ----

# Mantel test function
source(here("R Scripts", "Mantel Testing", "mantel_test.R"))

# Function to select donor from count matrix using donor label
select_donor <- function(df, donor_label) {
  
  donor_df <- df %>%
    select(gene_symbol, contains(donor_label), -contains("S18")) 
  colnames(donor_df) = gsub("_", "", colnames(donor_df))
  colnames(donor_df) = gsub(donor_label, "", colnames(donor_df))
  donor_df %<>% rename(S1 = S01) %>% rename(gene_symbol = genesymbol)
  
  return(donor_df)
}
# Function to create layer representations from donor-filtered df
create_layer_representations <- function(donor_df) {
  
  layer_df <- donor_df
  
  layer_df %<>%
    mutate(S4_weighted = S4 * 0.5) %>%
    mutate(S9_weighted = S9 * 0.5)
  
  return(layer_df)
}

# Create averaged He et al. tibble from the scaled values
create_layer_sum <- function(layer_df) {
  
  # Function to create a weighted average column
  wt_sum_col_fn <- function(layer_df, cols_to_sum) {
    df <- layer_df
    # Create wt.mean column - averages across columns
    sum_col <- df %>%
      select(all_of(cols_to_sum)) %>%
      data.table::transpose() %>%
      colSums()
    # Remove name attributes
    names(sum_col) <- NULL
    return(sum_col)
  }
  
  # Create tibble
  layer_sum_df <- tibble(
    gene_symbol = layer_df$gene_symbol,
    L1 = layer_df$S1,
    L2 = wt_sum_col_fn(layer_df, c("S2", "S3", "S4_weighted")),
    L3 = wt_sum_col_fn(layer_df, c("S4_weighted", "S5", "S6")),
    L4 = wt_sum_col_fn(layer_df, c("S7", "S8", "S9_weighted")),
    L5 = wt_sum_col_fn(layer_df, c("S9_weighted", "S10", "S11", "S12")),
    L6 = wt_sum_col_fn(layer_df, c("S13","S14", "S15", "S16")),
    WM = layer_df$S17
  )
  
  return(layer_sum_df)
}

# CPM normalize layer-wise tibble
CPM_normalize <- function(layer_sum_df, filter_for_low) {
  
  if (filter_for_low == TRUE) {
    
    CPM_df <- layer_sum_df %>%
      # Add one to counts to avoid taking cpm of 0
      mutate_at(c("L1", "L2", "L3", "L4", "L5", "L6", "WM"), ~. +1) %>%
      column_to_rownames(var = "gene_symbol") %>%
      cpm() %>%
      as.data.frame() %>%
      # Retain gene symbols, as filter removes rownames
      rownames_to_column(var = "gene_symbol") %>%
      # Filter out samples of CPM < 0.1
      filter_at(vars(-gene_symbol), all_vars(. > .1)) %>%
      column_to_rownames(var = "gene_symbol")
    names <- rownames(CPM_df)
    CPM_df %<>%
      # Take log2 of CPM
      map_df(log2) %>%
      select(L1, L2, L3, L4, L5, L6, WM) %>%
      # Take z-score (for app)
      t() %>% scale() %>% t() %>%
      as.data.frame() %>%
      add_column(gene_symbol = names) %>%
      select(gene_symbol, everything())
    
  } else {
    
    CPM_df <- layer_sum_df %>%
      # Remove gene_symbol column for cpm()
      select(-gene_symbol) %>% 
      # Add +1 to remove 0's for log transformation
      mutate(across(where(is.numeric), ~. +1)) %>%
      # CPM normalize with log = T
      cpm(log = T) %>% as.data.frame() %>%
      # Add back in gene symbols
      add_column(gene_symbol = layer_sum_df$gene_symbol) %>%
      select(gene_symbol, everything())
  }
  
  return(CPM_df)
}

## Data processing ----

# Read in He raw count matrix
He_count_matrix <- fread(file = here("Data", "raw_data", "He et al", "he_human_102_counts_matrix.csv"))

# Load Zeng layer markers
load(here("Data", "processed_data", "Zeng_dataset_long.Rdata"))
Zeng_marker_genes <- unique(Zeng_dataset_long$gene_symbol)

# Clean up column names
names(He_count_matrix) <- gsub(x = names(He_count_matrix), pattern = "\\-", replacement = "\\_") 

# Duplicate gene symbols exist - take the average across each section
He_count_matrix %<>%
  group_by(gene_symbol) %>%
  summarise(across(.fns = mean))

# Separate by DS1
He_DS1_matrix <- He_count_matrix %>%
  select(gene_symbol, contains("DS1")) %>%
  column_to_rownames(var = "gene_symbol") %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var = "sample_metadata")

# Clean up sample metadata column
He_DS1_matrix %<>%
  mutate(cleaned_metadata = str_replace(sample_metadata, ".*1_", ""),
         cleaned_metadata = str_replace(cleaned_metadata, "S1", "S01"),
         cleaned_metadata = str_replace(cleaned_metadata, "S010", "S10"),
         cleaned_metadata = str_replace(cleaned_metadata, "S011", "S11"),
         cleaned_metadata = str_replace(cleaned_metadata, "S012", "S12"),
         cleaned_metadata = str_replace(cleaned_metadata, "S013", "S13"),
         cleaned_metadata = str_replace(cleaned_metadata, "S014", "S14"),
         cleaned_metadata = str_replace(cleaned_metadata, "S015", "S15"),
         cleaned_metadata = str_replace(cleaned_metadata, "S016", "S16"),
         cleaned_metadata = str_replace(cleaned_metadata, "S017", "S17"),
         cleaned_metadata = str_replace(cleaned_metadata, "S018", "S18")) %>%
  select(-sample_metadata) %>% select(cleaned_metadata, everything())

# Transpose back to original dataframe (sample [col] by gene [row])
He_DS1_matrix %<>% 
  data.table::transpose(
    keep.names = "gene_symbol",
    make.names = names(He_DS1_matrix[1])
  ) %>% as.data.frame()

# Create donor-separated matrices
He_Human1 <- He_DS1_matrix %>%
  select("gene_symbol", "S01", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10",
         "S11", "S12", "S13", "S14", "S15", "S16", "S17") %>% 
  rename(S1 = S01)
He_Human1 %<>% 
  create_layer_representations() %>%
  create_layer_sum()

He_Human2 <- select_donor(He_DS1_matrix, "Human2") %>%
  create_layer_representations() %>%
  create_layer_sum()

He_Human3 <- select_donor(He_DS1_matrix, "Human3") %>%
  create_layer_representations() %>%
  create_layer_sum() 

He_Human4 <- select_donor(He_DS1_matrix, "Human4") %>%
  create_layer_representations() %>%
  create_layer_sum() 


# Create non-filtered CPM normalized data per donor
He_Human1_logCPM <- He_Human1 %>% CPM_normalize(filter_for_low = FALSE)
He_Human2_logCPM <- He_Human2 %>% CPM_normalize(filter_for_low = FALSE)
He_Human3_logCPM <- He_Human3 %>% CPM_normalize(filter_for_low = FALSE)
He_Human4_logCPM <- He_Human4 %>% CPM_normalize(filter_for_low = FALSE)

# Create filtered CPM normalized data per donor
He_Human1_logCPM_filtered <- He_Human1 %>% CPM_normalize(filter_for_low = TRUE)
He_Human2_logCPM_filtered <- He_Human2 %>% CPM_normalize(filter_for_low = TRUE)
He_Human3_logCPM_filtered <- He_Human3 %>% CPM_normalize(filter_for_low = TRUE)
He_Human4_logCPM_filtered <- He_Human4 %>% CPM_normalize(filter_for_low = TRUE)

# Create Zeng-filtered CPM normalized data per donor
He_Human1_Zeng_logCPM <- He_Human1_logCPM %>%
  filter(gene_symbol %in% Zeng_dataset_long$gene_symbol)
He_Human2_Zeng_logCPM <- He_Human2_logCPM %>%
  filter(gene_symbol %in% Zeng_dataset_long$gene_symbol)
He_Human3_Zeng_logCPM <- He_Human3_logCPM %>%
  filter(gene_symbol %in% Zeng_dataset_long$gene_symbol)
He_Human4_Zeng_logCPM <- He_Human4_logCPM %>%
  filter(gene_symbol %in% Zeng_dataset_long$gene_symbol)

# Clean up remaining DS1 data
rm(He_count_matrix, He_DS1_matrix, He_Human1, He_Human2, He_Human3, He_Human4)


### Mantel Testing ----

## Human 1 vs. Human 2
H1_H2 <- list()
for (type in c("logCPM", "logCPM_filtered", "Zeng_subset")) {
  H1_H2[type] <- mantel_test(He_Human1_logCPM, He_Human2_logCPM)
  H1_H2[type] <- mantel_test(He_Human1_logCPM, He_Human2_logCPM)
  H1_H2[type] <- mantel_test(He_Human1_logCPM, He_Human2_logCPM)
}
# logCPM
H1_H2[logCPM] <-  mantel_test(He_Human1_logCPM, He_Human2_logCPM)
# logCPM filtered for CPM > 0.1
mantel_test(He_Human_1_logCPM_filtered, He_Human2_logCPM_filtered)
# logCPM filtered for Zeng
H1_H2[['Zeng_subset']]<-  mantel_test(He_Human1_Zeng_logCPM, He_Human2_Zeng_logCPM)


## Human 1 vs. Human 3
# logCPM
mantel_test(He_Human1_logCPM, He_Human2_logCPM)
# logCPM filtered for CPM > 0.1
mantel_test(He_Human_1_logCPM_filtered, He_Human2_logCPM_filtered)
# logCPM filtered for Zeng
mantel_test(He_Human1_Zeng_logCPM, He_Human2_Zeng_logCPM)

## Human 1 vs. Human 4
# logCPM
mantel_test(He_Human1_logCPM, He_Human2_logCPM)
# logCPM filtered for CPM > 0.1
mantel_test(He_Human_1_logCPM_filtered, He_Human2_logCPM_filtered)
# logCPM filtered for Zeng
mantel_test(He_Human1_Zeng_logCPM, He_Human2_Zeng_logCPM)


## Human 2 vs. Human 3
# logCPM
mantel_test(He_Human1_logCPM, He_Human2_logCPM)
# logCPM filtered for CPM > 0.1
mantel_test(He_Human_1_logCPM_filtered, He_Human2_logCPM_filtered)
# logCPM filtered for Zeng
mantel_test(He_Human1_Zeng_logCPM, He_Human2_Zeng_logCPM)


## Human 2 vs. Human 4
# logCPM
mantel_test(He_Human1_logCPM, He_Human2_logCPM)
# logCPM filtered for CPM > 0.1
mantel_test(He_Human_1_logCPM_filtered, He_Human2_logCPM_filtered)
# logCPM filtered for Zeng
mantel_test(He_Human1_Zeng_logCPM, He_Human2_Zeng_logCPM)



## Human 3 vs. Human 4
# logCPM
mantel_test(He_Human1_logCPM, He_Human2_logCPM)
# logCPM filtered for CPM > 0.1
mantel_test(He_Human_1_logCPM_filtered, He_Human2_logCPM_filtered)
# logCPM filtered for Zeng
mantel_test(He_Human1_Zeng_logCPM, He_Human2_Zeng_logCPM)





