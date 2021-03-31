## Mantel tests for between scRNA-seq data ##

# Load libraries
library(dplyr)
library(magrittr)
library(conflicted) # Easily manage conflicting libraries
library(vegan) # Package for mantel test
library(here)
library(parallel) # Parallelize Mantel and WGCNA corr()
library(WGCNA) # Faster corr() than base

# Set conflicts
conflict_prefer('intersect', 'dplyr')
conflict_prefer("filter", "dplyr")
conflict_prefer("slice", "dplyr")

# Load AIBS_logCPM data
load(here("Data", "processed_data", "Allen_logCPM_dataset.Rdata"))

# Source mantel test function
source(here("R Scripts", "Mantel testing", "mantel_test.R"))

# Function to separate by cell type
separate_by_type <- function(AIBS_logCPM_data, cell_type) {
  
  df <- AIBS_logCPM_data
  
  cell_type_df <- df %>%
    filter(class_label == cell_type) %>%
    select(-class_label)
  
  return(cell_type_df)
}

# Create list of tibbles of all regions in layer(col) by gene(row) format
AIBS_cell_type <- list()
for (cell_type in unique(Allen_logCPM_dataset$class_label)) {
  AIBS_cell_type[[cell_type]] <- separate_by_type(Allen_logCPM_dataset, cell_type)
}


# Mantel tests between cell types ----

# GABAergic vs. Glutamatergic
# Mantel r: -0.1077, p < 0.001, n = 30,744
mantel_test(AIBS_cell_type$GABAergic, AIBS_cell_type$Glutamatergic)

# He logCPM vs Allen logCPM - Glutamatergic
# Mantel r: -0.1295, p < 0.001, n = 30,744
mantel_test(AIBS_cell_type$GABAergic, AIBS_cell_type$Glutamatergic)

# He logCPM vs Allen logCPM - Non-neuronal
# Mantel r: 0.09004, p < 0.001, n = 30,744
mantel_test(AIBS_cell_type$GABAergic, AIBS_cell_type$Glutamatergic)
