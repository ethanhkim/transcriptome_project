## Script to create randomly downsampled data from the AIBS snRNA-seq data.

# Load required libraries
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(purrr)
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
conflict_prefer("cpm", "edgeR")

# Read in data
MTG_cell_type_sum_count <- fread(here("Data", "raw_data", "Allen", "MTG_cell_type_sum_count.csv")) %>%
  select(-V1)
MTG_gene_sum_count <- fread(here("Data", "raw_data", "Allen", "MTG_gene_sum_count.csv")) %>%
  select(-V1)

GABA <- MTG_cell_type_sum_count %>%
  filter(class_label == 'GABAergic')
GLUT <- MTG_cell_type_sum_count %>%
  filter(class_label == 'Glutamatergic')
NONN <- MTG_cell_type_sum_count %>%
  filter(class_label == 'Non-neuronal')

GABA_layer_nuc_count <- GABA %>%
  group_by(cortical_layer_label) %>%
  count()
GLUT_layer_nuc_count <- GLUT %>%
  group_by(cortical_layer_label) %>%
  count()
NONN_layer_nuc_count <- NONN %>%
  group_by(cortical_layer_label) %>%
  count()

cell_type_nuc_count <- tibble(
  GABAergic = GABA_layer_nuc_count %>% pull(),
  Glutamatergic = GLUT_layer_nuc_count %>% pull(),
  Non_neuronal = NONN_layer_nuc_count %>% pull()
)

sample_nuclei <- function(df, layer_label, n_sample) {
  
  df %<>% filter(cortical_layer_label == layer_label)
  sampled_layer <- df[sample(nrow(df), n_sample),]
  return(sampled_layer)
}


GABA_sampled <- list()
GLUT_sampled <- list()
NONN_sampled <- list()
for (layer in c("L1", "L2", "L3", "L4", "L5", "L6")) {
  GABA_sampled[[layer]] <- sample_nuclei(GABA, layer, 381)
  GLUT_sampled[[layer]] <- sample_nuclei(GLUT, layer, 274)
  NONN_sampled[[layer]] <- sample_nuclei(NONN, layer, 125)
}

GABA_sampled_df <- rbindlist(GABA_sampled)
GLUT_sampled_df <- rbindlist(GLUT_sampled)
NONN_sampled_df <- rbindlist(NONN_sampled)

Allen_downsampledsampled_df <- rbind(GABA_sampled_df, GLUT_sampled_df, NONN_sampled_df)
Allen_downsampled_cell_sum_count <- sum_gene_count(Allen_downsampled_df, "cell_type")


write.csv(Allen_downsampled_cell_sum_count, 
          file = here("Data", "raw_data", "Allen", "Allen_downsampled_cell_sum_count.csv"))


Allen_downsampled_cell_count <- fread(here("Data", "raw_data", "Allen", "Allen_downsampled_cell_count.csv")) %>%
  select(-V1)

Allen_downsampled_logCPM_dataset <- Allen_downsampled_cell_count %>%
  # Add one to counts to avoid taking cpm of 0
  mutate_at(c("L1", "L2", "L3", "L4", "L5", "L6"), ~. +1) %>%
  unite(gene_class, c("gene_symbol", "class_label")) %>%
  column_to_rownames(var = "gene_class") %>%
  cpm(log = T) %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_class") %>%
  separate(gene_class, into = c("gene_symbol", "class_label"),
           sep = "_") %>%
  add_column(WM = NA) %>%
  select(gene_symbol, class_label, L1, L2, L3, L4, L5, L6, WM)


# Filtered data: CPM > 0.1
Allen_downsampled_logCPM_filtered_dataset <- Allen_downsampled_cell_count %>%
  # Add one to counts to avoid taking cpm of 0
  mutate_at(c("L1", "L2", "L3", "L4", "L5", "L6"), ~. +1) %>%
  unite(gene_class, c("gene_symbol", "class_label")) %>%
  column_to_rownames(var = "gene_class") %>%
  cpm() %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_class") %>%
  separate(gene_class, into = c("gene_symbol", "class_label"),
           sep = "_") %>%
  # Filter out samples of CPM < 0.1
  filter_at(vars(-gene_symbol, -class_label), all_vars(. > .1)) %>%
  unite(gene_class, c("gene_symbol", "class_label")) %>%
  column_to_rownames(var = "gene_class")
names <- rownames(Allen_downsampled_logCPM_filtered_dataset)
Allen_downsampled_logCPM_filtered_dataset %<>%
  # Take log2 of CPM
  map_df(log2) %>%
  select(L1, L2, L3, L4, L5, L6) %>%
  # Take z-score (for app)
  t() %>% scale() %>% t() %>%
  as.data.frame() %>%
  add_column(gene_class = names, WM = NA) %>%
  separate(gene_class, into = c("gene_symbol", "class_label"),
           sep = "_") %>%
  select(gene_symbol, class_label, L1, L2, L3, L4, L5, L6, WM)


save(Allen_downsampled_logCPM_dataset, file = here("Data", "processed_data", "Allen_downsampled_logCPM_dataset.Rdata"))
save(Allen_downsampled_logCPM_filtered_dataset, 
     file = here("Data", "processed_data", "Allen_downsampled_logCPM_filtered_dataset.Rdata"))