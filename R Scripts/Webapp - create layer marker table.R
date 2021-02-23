# Script to create laminar marker lookup table

library(tidyverse)
library(here)
library(magrittr)
library(spatialLIBD)
library(limma)
library(purrr)
library(conflicted)
library(readxl)

# Set conflicts
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("slice", "dplyr")

## He Layer markers ##
He_LM_Path <- here("Data", "raw_data", "He et al", "Supplementary Table 2.xlsx")
He_layer_markers <- read_xlsx(path = He_LM_Path) %>%
  select("Gene symbol", "Layer marker in human") %>%
  mutate_at(vars("Layer marker in human"), na_if, "No") %>%
  mutate_at(vars("Layer marker in human"), na_if, "NA") %>%
  rename(gene_symbol = "Gene symbol",  layer_marker = "Layer marker in human") %>%
  mutate(layer_marker = str_replace(layer_marker, "L", "")) %>%
  distinct(gene_symbol, .keep_all = T)

## Maynard layer markers ##
modeling_results <- fetch_data(type = "modeling_results")
Maynard_layer_markers <- as_tibble(modeling_results$enrichment) %>%
  distinct(gene, .keep_all = TRUE) %>%
  dplyr::rename(tstat_WM = "t_stat_WM", tstat_Layer1 = "t_stat_Layer1", tstat_Layer2 = "t_stat_Layer2", tstat_Layer3 = "t_stat_Layer3",
                tstat_Layer4 = "t_stat_Layer4", tstat_Layer5 = "t_stat_Layer5", tstat_Layer6 = "t_stat_Layer6", gene_symbol = "gene") %>%
  dplyr::select(-"ensembl", -starts_with("p")) %>%
  filter(gene_symbol %in% modeling_results$enrichment$gene) %>%
  pivot_longer(
    cols = -gene_symbol,
    names_to = c(".value", "initial_layer_marker"),
    names_sep = "_"
  ) %>%
  group_by(gene_symbol) %>% 
  slice(which.max(tstat)) %>%
  mutate(fdr_true_false = ifelse(test = (fdr < 0.1),
                                 yes = "yes",
                                 no = "no")) %>%
  mutate(layer_marker = ifelse(test = (fdr_true_false == "yes"),
                               yes = initial_layer_marker,
                               no = NA)) %>%
  select(gene_symbol, layer_marker) %>%
  mutate(layer_marker = str_replace(layer_marker, "Layer", "")) %>%
  mutate(layer_marker = str_replace(layer_marker, "WM", "7")) %>%
  distinct(gene_symbol, .keep_all = T)

# Create lookup table for all layer markers
layer_marker_table <- left_join(Maynard_layer_markers, He_layer_markers, by = "gene_symbol") %>%
  rename(He = layer_marker.x) %>%
  rename(Maynard = layer_marker.y)

# Export lookup table
save(layer_marker_table, file = here("Data", "processed_data", "layer_marker_table.Rdata"))
