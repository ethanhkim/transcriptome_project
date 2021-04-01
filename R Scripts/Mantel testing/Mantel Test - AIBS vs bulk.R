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

# He logCPM vs Allen logCPM - GABAergic
# Mantel r: -0.1077, p < 0.001, n = 30,744
He_AIBS_logCPM <- list()
for (type in c("GABAergic", "Glutamatergic", "Non-neuronal")) {
  He_AIBS_logCPM[[type]] <- mantel_test(He_DS1_logCPM_dataset,
                                        Allen_cell_type_df[[type]])
}

# He logCPM vs Allen logCPM - Glutamatergic
# Mantel r: -0.1295, p < 0.001, n = 30,744
He_AIBS_GLUT <- mantel_test(He_DS1_logCPM_dataset, Allen_cell_type_df$GABAergic)

# He logCPM vs Allen logCPM - Non-neuronal
# Mantel r: 0.09004, p < 0.001, n = 30,744
He_AIBS_NONN <- mantel_test(He_DS1_logCPM_dataset, Allen_NONN)

## He et al vs. Allen - filtered for CPM > 0.1 ## ----

# He  vs Allen - GABAergic
# Mantel r: -0.006967, p < 0.001, n = 13,114
mantel_test(He_DS1_logCPM_filtered_dataset, 
            Allen_GABA_filtered)

# He  vs Allen - Glutamatergic
# Mantel r: -0.005718, p < 0.001, n = 13,114
mantel_test(He_DS1_logCPM_filtered_dataset, 
            Allen_GLUT_filtered)

# He  vs Allen - Non-neuronal
# Mantel r: 0.05965, p < 0.001, n = 11,645
mantel_test(He_DS1_logCPM_filtered_dataset, 
            Allen_NONN_filtered)

## He et al vs. Allen - filtered for Zeng ----

# He vs. Allen - GABAergic
# Mantel r: 0.03838, p < 0.001, n = 975
mantel_test(He_Zeng, Allen_GABA_Zeng)

# He vs. Allen - Glutamatergic
# Mantel r: 0.03376, p < 0.001, n = 975
mantel_test(He_Zeng, Allen_GLUT_Zeng)

# He vs. Allen - Non-neuronal
# Mantel r: 0.07922, p < 0.001, n = 975
mantel_test(He_Zeng, Allen_NONN_Zeng)


## Maynard et al vs. Allen ## ----

# Maynard logCPM vs Allen logCPM - GABAergic
# Mantel r: -0.01428, p < 0.001, n = 17,197
mantel_test(Maynard_logCPM_dataset, Allen_GABA)

# Maynard logCPM vs Allen logCPM - Glutamatergic
# Mantel r: -0.0143, p < 0.001, n = 17,197
mantel_test(Maynard_logCPM_dataset, Allen_GLUT)

# Maynard logCPM vs Allen logCPM - Non-neuronal
# Mantel r: 0.02628, p < 0.001, n = 17,197
mantel_test(Maynard_logCPM_dataset, Allen_NONN)


## Maynard et al vs. Allen - filtered for CPM > 0.1 ## ----

# Maynard  vs Allen - GABAergic
# Mantel r: 0.002822, p < 0.001, n = 12,913
mantel_test(Maynard_logCPM_filtered_dataset, 
            Allen_GABA_filtered)

# Maynard  vs Allen - Glutamatergic
# Mantel r: 0.0097692, p < 0.001, n = 13,614
mantel_test(Maynard_logCPM_filtered_dataset, 
            Allen_GLUT_filtered)

# Maynard  vs Allen - Non-neuronal
# Mantel r: 0.04014, p < 0.001, n = 11,471
mantel_test(Maynard_logCPM_filtered_dataset, 
            Allen_NONN_filtered)

## Maynard et al vs. Allen - filtered for Zeng ## ----

# Maynard vs. Allen - GABAergic
# Mantel r: 0.03121, p < 0.001, n = 929
mantel_test(Maynard_Zeng, Allen_GABA_Zeng)

# Maynard vs. Allen - Glutamatergic
# Mantel r: 0.064, p < 0.001, n = 929
mantel_test(Maynard_Zeng, Allen_GLUT_Zeng)

# Maynard vs. Allen - Non-neuronal
# Mantel r: 0.06274, p < 0.001, n = 929
Maynard_AIBS_Zeng_NONN <- mantel_test(Maynard_Zeng, Allen_NONN_Zeng)

mantel_test(Maynard_logCPM_dataset, He_DS1_logCPM_dataset)

