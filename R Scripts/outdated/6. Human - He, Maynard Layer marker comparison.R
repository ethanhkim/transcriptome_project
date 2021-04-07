library(tidyverse)
library(here)
library(magrittr)
library(spatialLIBD)
library(limma)
library(purrr)
library(dplyr)

## He Layer markers ##
He_LM_Path <- here("Data", "He et al", "Supplementary Table 2.xlsx")
He_Layer_markers <- read_xlsx(path = He_LM_Path) %>%
  dplyr::select("Gene symbol", "Layer marker in human") %>%
  mutate_at(vars("Layer marker in human"), na_if, "No") %>%
  mutate_at(vars("Layer marker in human"), na_if, "NA") %>%
  dplyr::rename(gene_symbol = "Gene symbol",  layer_marker = "Layer marker in human")

## Maynard layer markers ##
modeling_results <- fetch_data(type = "modeling_results")
sce_layer_data <- sce_layer@rowRanges@elementMetadata$gene_name
Maynard_modeling <- as_tibble(modeling_results$enrichment) %>%
  distinct(gene, .keep_all = TRUE) %>%
  dplyr::rename(tstat_WM = "t_stat_WM", tstat_Layer1 = "t_stat_Layer1", tstat_Layer2 = "t_stat_Layer2", tstat_Layer3 = "t_stat_Layer3",
                tstat_Layer4 = "t_stat_Layer4", tstat_Layer5 = "t_stat_Layer5", tstat_Layer6 = "t_stat_Layer6", gene_symbol = "gene") %>%
  dplyr::select(-"ensembl", -starts_with("p")) %>%
  filter(gene_symbol %in% sce_layer_data)

## Create common gene list between He, Maynard
Compared_gene_list <- intersect(He_Layer_markers$gene_symbol, Maynard_layer_enrichment$gene_symbol)
He_Layer_Markers <- as_tibble(He_Layer_markers[He_Layer_markers$gene_symbol %in% Compared_gene_list,]) %>%
  arrange(factor(gene_symbol, levels = Compared_gene_list)) %>%
  distinct(gene_symbol, .keep_all = TRUE)
Maynard_Layer_Markers <- Maynard_layer_enrichment[Maynard_layer_enrichment$gene_symbol %in% Compared_gene_list,]

## Create tibble of layer markers ##
Compared_Layer_Markers <- tibble(
  gene_symbol = Compared_gene_list,
  He_Layer_Markers = He_Layer_Markers$layer_marker,
  Maynard_Layer_Markers = Maynard_Layer_Markers$layer_marker) 

## Re-format ##
Compared_Layer_Markers %<>%
  mutate(He_Layer_Markers = gsub("L", "", He_Layer_Markers)) %>%
  mutate(Maynard_Layer_Markers = gsub("Layer", "", Maynard_Layer_Markers)) %>%
  mutate(Maynard_Layer_Markers = gsub("WM", "7", Maynard_Layer_Markers))

Compared_Layer_Markers_filtered <- Compared_Layer_Markers %>%
  na.exclude()
Compared_Layer_Markers_filtered$He_Layer_Markers <- as.numeric(as.character(Compared_Layer_Markers_filtered$He_Layer_Markers))
Compared_Layer_Markers_filtered$Maynard_Layer_Markers <- as.numeric(as.character(Compared_Layer_Markers_filtered$Maynard_Layer_Markers))

Compared_Layer_Markers_filtered %<>%
  mutate(Maynard_Layer_Markers = ifelse(
    test = Maynard_Layer_Markers == 7,
    yes = 6,
    no = Maynard_Layer_Markers
  )) %<>%
  mutate(match = ifelse(
    test = (abs(Maynard_Layer_Markers - He_Layer_Markers) < 2),
    yes = "YES",
    no = "NO"
  )) 

Compared_Layer_Markers_filtered %>%
  group_by(match) %>%
  summarize(n())

Compared_Layer_Markers_filtered %>%
  group_by(Maynard_Layer_Markers) %>%
  summarize(n())
