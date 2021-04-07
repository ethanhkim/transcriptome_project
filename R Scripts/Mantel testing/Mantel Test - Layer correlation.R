## Mantel Test - Layer-to-layer comparisons ##

# Load libraries
library(dplyr)
library(magrittr)
library(tibble)
library(conflicted) # Easily manage conflicting libraries
library(vegan) # Package for mantel test
library(here)
library(parallel) # Parallelize Mantel and WGCNA corr()
library(WGCNA) # Faster corr() than base
library(gplots)

# Set conflicts
conflict_prefer('intersect', 'dplyr')
conflict_prefer("filter", "dplyr")
conflict_prefer("slice", "dplyr")

# Source mantel test function
source(here("R Scripts", "Mantel testing", "mantel_test.R"))

# Load in data
load(here("Data", "processed_data", "He_DS1_logCPM_dataset.Rdata"))
load(here("Data", "processed_data", "He_DS1_logCPM_filtered_dataset.Rdata"))
load(here("Data", "processed_data", "Maynard_logCPM_dataset.Rdata"))
load(here("Data", "processed_data", "Maynard_logCPM_filtered_dataset.Rdata"))
load(here("Data", "processed_data", "Allen_logCPM_dataset.Rdata"))
load(here("Data", "processed_data", "Allen_logCPM_filtered_dataset.Rdata"))
load(here("Data", "processed_data", "Allen_gene_logCPM_dataset.Rdata"))
load(here("Data", "processed_data", "Allen_gene_logCPM_filtered_dataset.Rdata"))

# Format data

# Create correlation matrices for each datasets
correlation_matrices <- list()  

# Create heatmaps for each of the layer-to-layer correlations
