#### Allen vs. Maynard vs. He ####

# Load libraries
library(dplyr)
library(magrittr)
library(tibble)
library(conflicted) # Easily manage conflicting libraries
library(vegan) # Package for mantel test
library(here)
library(parallel) # Parallelize Mantel and WGCNA corr()
library(WGCNA) # Faster corr() than base

# Set conflicts
conflict_prefer('intersect', 'dplyr')
conflict_prefer("filter", "dplyr")
conflict_prefer("slice", "dplyr")

# Source mantel test function
source(here("R Scripts", "Mantel testing", "mantel_test.R"))

# Load in data
load(here("Data", "processed_data", "He_DS1_logCPM_dataset.Rdata"))
load(here("Data", "processed_data", "Maynard_logCPM_dataset.Rdata"))
load(here("Data", "processed_data", "He_DS1_logCPM_filtered_dataset.Rdata"))
load(here("Data", "processed_data", "Maynard_logCPM_filtered_dataset.Rdata"))
load(here("Data", "processed_data", "Allen_logCPM_dataset.Rdata"))
load(here("Data", "processed_data", "Allen_logCPM_filtered_dataset.Rdata"))

# Zeng markers
load(here("Data", "processed_data", "Zeng_dataset_long.Rdata"))
Zeng_marker_genes <- unique(Zeng_dataset_long$gene_symbol)

mantel_test(He_DS1_logCPM_dataset, Maynard_logCPM_dataset)

# Remove WM from He and Maynard
He_DS1_logCPM_dataset %<>% select(-WM)
He_DS1_logCPM_filtered_dataset %<>% select(-WM)
Maynard_logCPM_dataset %<>% select(-WM)
Maynard_logCPM_filtered_dataset %<>% select(-WM)

# Subset He and Maynard by Zeng markers
He_Zeng <- He_DS1_logCPM_dataset %>%
  filter(gene_symbol %in% Zeng_marker_genes)
Maynard_Zeng <- Maynard_logCPM_dataset %>%
  filter(gene_symbol %in% Zeng_marker_genes)

# Cell-type specific AIBS
Allen_cell_type_df <- list()
for (type in c("GABAergic", "Glutamatergic", "Non-neuronal")) {
  Allen_cell_type_df[[type]] <- Allen_logCPM_dataset %>%
    filter(class_label == type) %>%
    select(-class_label, -WM)
}

# Cell-type specific AIBS - filtered for CPM > 0.1
Allen_cell_type_filtered_df <- list()
for (type in c("GABAergic", "Glutamatergic", "Non-neuronal")) {
  Allen_cell_type_filtered_df[[type]] <- Allen_logCPM_filtered_dataset %>%
    filter(class_label == type) %>%
    select(-class_label, -WM)
}

# Cell-type specific AIBS - filtered for Zeng
Allen_cell_type_Zeng <- list()
for (type in c("GABAergic", "Glutamatergic", "Non-neuronal")) {
  Allen_cell_type_Zeng[[type]] <- Allen_logCPM_dataset %>%
    filter(class_label == type) %>%
    filter(gene_symbol %in% Zeng_marker_genes) %>%
    select(-class_label, -WM)
}



# Mantel tests b/w He, Maynard, Allen ----

## He et al vs. Allen ## ---- 

## He logCPM vs Allen logCPM ----
# GABAergic -     Mantel r: -0.1194, p < 0.001, n = 30,744
# Glutamatergic - Mantel r: -0.1526, p < 0.001, n = 30,744
# Non-neuronal -  Mantel r: 0.1025, p < 0.001, n = 30,744
He_AIBS_logCPM <- list()
for (type in c("GABAergic", "Glutamatergic", "Non-neuronal")) {
  He_AIBS_logCPM[[type]] <- mantel_test(He_DS1_logCPM_dataset,
                                        Allen_cell_type_df[[type]])
}
for (type in c("GABAergic", "Glutamatergic", "Non-neuronal")) {
  print(He_AIBS_logCPM[[type]]$statistic)
}

## He logCPM vs Allen logCPM - filtered for CPM > 0.1  ---- 
# GABAergic -     Mantel r: -0.005578, p < 0.001, n = 13,114
# Glutamatergic - Mantel r: 0.003160, p < 0.001, n = 13,714
# Non-neuronal -  Mantel r: 0.06069, p < 0.001, n = 11,645
He_AIBS_logCPM_filtered <- list()
for (type in c("GABAergic", "Glutamatergic", "Non-neuronal")) {
  He_AIBS_logCPM_filtered[[type]] <- mantel_test(He_DS1_logCPM_filtered_dataset,
                                        Allen_cell_type_filtered_df[[type]])
}
for (type in c("GABAergic", "Glutamatergic", "Non-neuronal")) {
  print(He_AIBS_logCPM_filtered[[type]]$statistic)
}

## He logCPM vs Allen logCPM - Zeng subset ----
# GABAergic -     Mantel r: 0.02548, p < 0.001, n = 13,114
# Glutamatergic - Mantel r: 0.03520, p < 0.001, n = 13,114
# Non-neuronal -  Mantel r: 0.07255, p < 0.001, n = 11,645
He_AIBS_Zeng <- list()
for (type in c("GABAergic", "Glutamatergic", "Non-neuronal")) {
  He_AIBS_Zeng[[type]] <- mantel_test(He_Zeng, Allen_cell_type_Zeng[[type]])
}
for (type in c("GABAergic", "Glutamatergic", "Non-neuronal")) {
  print(He_AIBS_Zeng[[type]]$statistic)
}


## Maynard et al vs. Allen ## ----

## Maynard logCPM vs Allen logCPM ----
# GABAergic -     Mantel r: -0.01683, p < 0.001, n = 17,197
# Glutamatergic - Mantel r: -0.01113, p < 0.001, n = 17,197
# Non-neuronal -  Mantel r: 0.02406, p < 0.001, n = 17,197
Maynard_AIBS_logCPM <- list()
for (type in c("GABAergic", "Glutamatergic", "Non-neuronal")) {
  Maynard_AIBS_logCPM[[type]] <- mantel_test(Maynard_logCPM_dataset,
                                        Allen_cell_type_df[[type]])
}
for (type in c("GABAergic", "Glutamatergic", "Non-neuronal")) {
  print(Maynard_AIBS_logCPM[[type]]$statistic)
}


## Maynard et al vs. Allen - filtered for CPM > 0.1 ## ----
# GABAergic -     Mantel r: 0.006462, p < 0.001, n = 12,913
# Glutamatergic - Mantel r: 0.002535, p < 0.001, n = 13,614
# Non-neuronal -  Mantel r: 0.03572, p < 0.001, n = 11,471
Maynard_AIBS_logCPM_filtered <- list()
for (type in c("GABAergic", "Glutamatergic", "Non-neuronal")) {
  Maynard_AIBS_logCPM_filtered[[type]] <- mantel_test(Maynard_logCPM_filtered_dataset,
                                             Allen_cell_type_filtered_df[[type]])
}
for (type in c("GABAergic", "Glutamatergic", "Non-neuronal")) {
  print(Maynard_AIBS_logCPM_filtered[[type]]$statistic)
}


## Maynard et al vs. Allen - filtered for Zeng ## ----
# GABAergic -     Mantel r: 0.008808, p < 0.001, n = 929
# Glutamatergic - Mantel r: 0.07633, p < 0.001, n = 929
# Non-neuronal -  Mantel r: 0.04006, p < 0.001, n = 929
Maynard_AIBS_Zeng <- list()
for (type in c("GABAergic", "Glutamatergic", "Non-neuronal")) {
  Maynard_AIBS_Zeng[[type]] <- mantel_test(Maynard_Zeng, Allen_cell_type_Zeng[[type]])
}
for (type in c("GABAergic", "Glutamatergic", "Non-neuronal")) {
  print(Maynard_AIBS_Zeng[[type]]$statistic)
}



