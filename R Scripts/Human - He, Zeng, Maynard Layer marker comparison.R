library(tidyverse)
library(here)
library(magrittr)
library(readxl)
library(here)
library(spatialLIBD)
library(limma)
library(purrr)
library(corrr)

## Zeng Layer markers ##
Zeng_Layer_markers <- Zeng_dataset %>%
  dplyr::select("gene_symbol", "marker_annotation") %>%
  mutate(marker_annotation = gsub("layer", "", marker_annotation)) %>%
  mutate(marker_annotation = gsub("interneuron", "", marker_annotation)) %>%
  mutate(marker_annotation = gsub("astrocyte", "", marker_annotation)) %>%
  mutate(marker_annotation = gsub("astrocyte?", "", marker_annotation)) %>%
  mutate(marker_annotation = gsub("laminar", "", marker_annotation)) %>%
  mutate(marker_annotation = gsub(" ", "", marker_annotation)) %>%
  mutate(marker_annotation = gsub("VEC", "", marker_annotation)) %>%
  mutate(marker_annotation = gsub("or", "", marker_annotation)) %>%
  mutate(marker_annotation = gsub("others", "", marker_annotation)) %>%
  mutate(marker_annotation = gsub("\\+", "", marker_annotation)) %>%
  mutate(marker_annotation = gsub("oligodendrocyte", "", marker_annotation)) %>%
  mutate(marker_annotation = gsub("5a", "5", marker_annotation)) %>%
  mutate(marker_annotation = gsub("6b", "6", marker_annotation)) %>%
  mutate(marker_annotation = gsub("4c", "4", marker_annotation)) %>%
  na_if("")

## He Layer markers ##
He_LM_Path <- here("Data", "He et al", "Supplementary Table 2.xlsx")
He_Layer_markers <- read_xlsx(path = He_LM_Path) %>%
  dplyr::select("Gene symbol", "Layer marker in human") %>%
  mutate_at(vars("Layer marker in human"), na_if, "No") %>%
  mutate_at(vars("Layer marker in human"), na_if, "NA") %>%
  dplyr::rename(gene_symbol = "Gene symbol",  layer_marker = "Layer marker in human")

## Maynard layer markers ##
modeling_results <- fetch_data(type = "modeling_results")
Maynard_modeling <- as_tibble(modeling_results$enrichment) %>%
  distinct(gene, .keep_all = TRUE) %>%
  dplyr::rename(tstat_WM = "t_stat_WM", tstat_Layer1 = "t_stat_Layer1", tstat_Layer2 = "t_stat_Layer2", tstat_Layer3 = "t_stat_Layer3",
                tstat_Layer4 = "t_stat_Layer4", tstat_Layer5 = "t_stat_Layer5", tstat_Layer6 = "t_stat_Layer6", gene_symbol = "gene") %>%
  dplyr::select(-"ensembl", -starts_with("p")) %>%
  filter(gene_symbol %in% Maynard_dataset_average$gene_symbol)

Maynard_layer_enrichment <- Maynard_modeling %>%
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
  mutate(layer_marker_label = ifelse(test = (fdr_true_false == "yes"),
                                     yes = "*",
                                     no = NA)) %>%
  mutate(layer_marker_label_1 = ifelse(test = (layer_marker == "Layer1"),
                                       yes = "Layer_1",
                                       no = NA)) %>%
  mutate(layer_marker_label_2 = ifelse(test = (layer_marker == "Layer2"),
                                       yes = "Layer_2",
                                       no = NA)) %>%
  mutate(layer_marker_label_3 = ifelse(test = (layer_marker == "Layer3"),
                                       yes = "Layer_3",
                                       no = NA)) %>%
  mutate(layer_marker_label_4 = ifelse(test = (layer_marker == "Layer4"),
                                       yes = "Layer_4",
                                       no = NA)) %>%
  mutate(layer_marker_label_5 = ifelse(test = (layer_marker == "Layer5"),
                                       yes = "Layer_5",
                                       no = NA)) %>%
  mutate(layer_marker_label_6 = ifelse(test = (layer_marker == "Layer6"),
                                       yes = "Layer_6",
                                       no = NA)) %>%
  mutate(layer_marker_label_WM = ifelse(test = (layer_marker == "WM"),
                                        yes = "WM",
                                        no = NA)) %>%
  dplyr::select(-tstat, -fdr)

## Create common gene list between He, Zeng, Maynard
Compared_gene_list <- intersect(intersect(Zeng_Layer_markers$gene_symbol, He_Layer_markers$gene_symbol), Maynard_layer_enrichment$gene_symbol)
He_Layer_Markers <- as_tibble(He_Layer_markers[He_Layer_markers$gene_symbol %in% Compared_gene_list,]) %>%
  arrange(factor(gene_symbol, levels = Compared_gene_list))
Zeng_Layer_Markers <- Zeng_Layer_markers[Zeng_Layer_markers$gene_symbol %in% Compared_gene_list,]
Maynard_Layer_Markers <- Maynard_layer_enrichment[Maynard_layer_enrichment$gene_symbol %in% Compared_gene_list,]


## Create tibble of layer markers ##
Compared_Layer_Markers <- tibble(
  gene_symbol = Compared_gene_list,
  He_Layer_Markers = He_Layer_Markers$layer_marker,
  Zeng_Layer_Markers = Zeng_Layer_Markers$layer_marker,
  Maynard_Layer_Markers = Maynard_Layer_Markers$layer_marker) 

## Re-format ##
Compared_Layer_Markers %<>%
  mutate(Zeng_Layer_Markers = gsub("\\?", NA, Zeng_Layer_Markers)) %>%
  mutate(Zeng_Layer_Markers = gsub("\\/1", "1", Zeng_Layer_Markers)) %>%
  mutate(Zeng_Layer_Markers = gsub("1\\/2\\/\\/", "1\\/2", Zeng_Layer_Markers)) %>%
  mutate(Zeng_Layer_Markers = gsub("\\/2\\/3", "2\\/3", Zeng_Layer_Markers)) %>%
  mutate(Zeng_Layer_Markers = gsub("34\\/5\\/6", "3\\/4\\/5\\/6", Zeng_Layer_Markers)) %>%
  mutate(Zeng_Layer_Markers = gsub("\\/2\\/3", "2\\/3", Zeng_Layer_Markers)) %>%
  mutate(Zeng_Layer_Markers = gsub("6\\/\\/", "6", Zeng_Layer_Markers)) %>%
  mutate(Zeng_Layer_Markers = gsub("6\\/6", "6", Zeng_Layer_Markers)) %>%
  mutate(Zeng_Layer_Markers = gsub("2\\/3\\/5\\/6", "2,3,5,6", Zeng_Layer_Markers)) %>%
  mutate(Zeng_Layer_Markers = gsub("3\\/4\\/5\\/6", "3,4,5,6", Zeng_Layer_Markers)) %>%
  mutate(Zeng_Layer_Markers = gsub("3\\/5\\/6", "3,5,6", Zeng_Layer_Markers)) %>%
  mutate(Zeng_Layer_Markers = gsub("2\\/3\\/6", "2,3,6", Zeng_Layer_Markers)) %>%
  mutate(Zeng_Layer_Markers = gsub("2\\/3\\/4", "2,3,4", Zeng_Layer_Markers)) %>%
  mutate(Zeng_Layer_Markers = gsub("2\\/3\\/5", "2,3,5", Zeng_Layer_Markers)) %>%
  mutate(Zeng_Layer_Markers = gsub("4\\/5\\/6", "3,5,6", Zeng_Layer_Markers)) %>%
  mutate(Zeng_Layer_Markers = gsub("1\\/2\\/6", "1,2,6", Zeng_Layer_Markers)) %>%
  mutate(Zeng_Layer_Markers = gsub("\\/5\\/6", "5,6", Zeng_Layer_Markers)) %>%
  mutate(Zeng_Layer_Markers = gsub("\\/4\\/5", "4,5", Zeng_Layer_Markers)) %>%
  mutate(Zeng_Layer_Markers = gsub("5\\/6", "5,6", Zeng_Layer_Markers)) %>%
  mutate(Zeng_Layer_Markers = gsub("3\\/5", "3,5", Zeng_Layer_Markers)) %>%
  mutate(Zeng_Layer_Markers = gsub("2\\/3", "2,3", Zeng_Layer_Markers)) %>%
  mutate(Zeng_Layer_Markers = gsub("1\\/2", "1,2", Zeng_Layer_Markers)) %>%
  mutate(Zeng_Layer_Markers = gsub("4\\/6", "4,6", Zeng_Layer_Markers)) %>%
  mutate(Zeng_Layer_Markers = gsub("3\\/4", "3,4", Zeng_Layer_Markers)) %>%
  mutate(Zeng_Layer_Markers = gsub("/2", "2", Zeng_Layer_Markers)) %>%
  mutate(Zeng_Layer_Markers = gsub("/4", "4", Zeng_Layer_Markers)) %>%
  mutate(Zeng_Layer_Markers = gsub("/6", "6", Zeng_Layer_Markers)) %>%
  mutate(Zeng_Layer_Markers = gsub("\\/", NA, Zeng_Layer_Markers)) %>%
  mutate(He_Layer_Markers = gsub("L", "", He_Layer_Markers)) %>%
  mutate(Maynard_Layer_Markers = gsub("Layer", "", Maynard_Layer_Markers)) %>%
  mutate(Maynard_Layer_Markers = gsub("WM", "7", Maynard_Layer_Markers)) %>%
  mutate(Zeng_Layer_Markers = strsplit(as.character(Zeng_Layer_Markers), "[,]")) %>% 
  unnest(Zeng_Layer_Markers)


## Goal: create columns of gene_symbol, layer, dataset

Compared_Layer_Markers_formatted <- Compared_Layer_Markers %>%
  rename(LayerMarkers_He = "He_Layer_Markers", LayerMarkers_Zeng = "Zeng_Layer_Markers", LayerMarkers_Maynard = "Maynard_Layer_Markers") %>%
  pivot_longer(
    cols = -gene_symbol,
    names_to = c(".value", "dataset"),
    names_sep = "_"
  ) %>%
  rename(layer_marker = "LayerMarkers", source_dataset = "dataset") %>%
  dplyr::select("gene_symbol", "layer_marker", "source_dataset") %>%
  mutate(layer_marker = strsplit(as.character(layer_marker), "[,]")) %>% 
  unnest(layer_marker) %>%
  as_tibble()

## Create common gene list between He, Maynard
Compared_gene_list_He_Maynard <- intersect(He_Layer_markers$gene_symbol, Maynard_layer_enrichment$gene_symbol)
He_Layer_Markers <- as_tibble(He_Layer_markers[He_Layer_markers$gene_symbol %in% Compared_gene_list,]) %>%
  arrange(factor(gene_symbol, levels = Compared_gene_list))
Maynard_Layer_Markers <- Maynard_layer_enrichment[Maynard_layer_enrichment$gene_symbol %in% Compared_gene_list,]



Compared_Layer_Markers_filtered <- Compared_Layer_Markers %>%
  na.exclude()
Compared_Layer_Markers_filtered$He_Layer_Markers <- as.numeric(as.character(Compared_Layer_Markers_filtered$He_Layer_Markers))
Compared_Layer_Markers_filtered$Zeng_Layer_Markers <- as.numeric(as.character(Compared_Layer_Markers_filtered$Zeng_Layer_Markers))
Compared_Layer_Markers_filtered$Maynard_Layer_Markers <- as.numeric(as.character(Compared_Layer_Markers_filtered$Maynard_Layer_Markers))

Compared_Layer_Markers_filtered %<>%
  mutate(match = ifelse(
    test = (
      (((Maynard_Layer_Markers == Zeng_Layer_Markers) & (He_Layer_Markers == Zeng_Layer_Markers))
      )),
    yes = "YES",
    no = "NO"
  )) %>%
  filter(str_detect(match, "YES"), .keep_all = TRUE)


write_csv(Compared_Layer_Markers, './R Scripts/export_data/Compared_Layer_Markers.csv')
save(Compared_Layer_Markers, file = './R Scripts/export_data/Compared_Layer_markers.Rdata')

