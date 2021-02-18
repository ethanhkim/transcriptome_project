## Extract layer markers ##

# Required libraries #
library(tidyverse)
library(spatialLIBD)
library(conflicted)
library(readxl)

conflict_prefer('slice', 'dplyr')
conflict_prefer('filter', 'dplyr')

## Maynard layer markers ##
modeling_results <- fetch_data(type = "modeling_results")
Maynard_modeling_results <- as_tibble(modeling_results$enrichment) %>%
  distinct(gene, .keep_all = TRUE) %>%
  dplyr::rename(tstat_WM = "t_stat_WM", tstat_Layer1 = "t_stat_Layer1", 
                tstat_Layer2 = "t_stat_Layer2", tstat_Layer3 = "t_stat_Layer3",
                tstat_Layer4 = "t_stat_Layer4", tstat_Layer5 = "t_stat_Layer5", 
                tstat_Layer6 = "t_stat_Layer6", pvalue_Layer1 = "p_value_Layer1",
                pvalue_Layer2 = "p_value_Layer2", pvalue_Layer3 = "p_value_Layer3",
                pvalue_Layer4 = "p_value_Layer4", pvalue_Layer5 = "p_value_Layer5",
                pvalue_Layer6 = "p_value_Layer6", pvalue_WM = "p_value_WM",
                gene_symbol = "gene") %>%
  dplyr::select(-"ensembl") %>%
  filter(gene_symbol %in% Maynard_logCPM_averaged$gene_symbol)

Maynard_layer_enrichment <- Maynard_modeling_results %>%
  pivot_longer(
    cols = -gene_symbol,
    names_to = c(".value", "cortical_layer"),
    names_sep = "_"
  ) %>%
  group_by(gene_symbol) %>% 
  slice(which.max(tstat)) %>%
  ungroup() %>%
  mutate(fdr_pass = ifelse(test = (fdr < 0.1),
                           yes = "pass",
                           no = "fail")) %>%
  mutate(pval_pass = ifelse(test = (pvalue < 0.05),
                            yes = "pass",
                            no = "fail")) %>%
  mutate(layer_marker = ifelse(test = (fdr_pass == "pass") & (pval_pass == "pass"),
                               yes = cortical_layer,
                               no = NA)) %>%
  dplyr::select(gene_symbol, layer_marker) 

## He et al Layer Markers ##

He_LM_Path <- here("Data", "raw_data", "He et al", "Supplementary Table 2.xlsx")
He_layer_marker <- read_xlsx(path = He_LM_Path) %>%
  dplyr::select("Gene symbol", "Layer marker in human") %>%
  mutate_at(vars("Layer marker in human"), na_if, "No") %>%
  mutate_at(vars("Layer marker in human"), na_if, "NA") %>%
  dplyr::rename(gene_symbol = "Gene symbol",  layer_marker = "Layer marker in human") %>%
  mutate(layer_marker = str_replace(layer_marker, "L", "Layer"))


## Zeng et al Layer Markers ##

Zeng_marker_annotation <- Zeng_dataset %>%
  dplyr::select("gene_symbol", "marker_annotation") %>%
  rename(layer_marker_Zeng = "marker_annotation") %>%
  mutate(layer_marker_Zeng = gsub("layer", "", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("interneuron", "", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("astrocyte", "", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("astrocyte?", "", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("laminar", "", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub(" ", "", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("VEC", "", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("or", "", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("others", "", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("\\+", "", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("oligodendrocyte", "", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("5a", "5", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("6b", "6", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("4c", "4", layer_marker_Zeng)) %>%
  na_if("") %>%
  mutate(layer_marker_Zeng = gsub("\\?", NA, layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("\\/1", "1", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("1\\/2\\/\\/", "1\\/2", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("\\/2\\/3", "2\\/3", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("34\\/5\\/6", "3\\/4\\/5\\/6", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("\\/2\\/3", "2\\/3", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("6\\/\\/", "6", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("6\\/6", "6", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("2\\/3\\/5\\/6", "2,3,5,6", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("3\\/4\\/5\\/6", "3,4,5,6", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("3\\/5\\/6", "3,5,6", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("2\\/3\\/6", "2,3,6", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("2\\/3\\/4", "2,3,4", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("2\\/3\\/5", "2,3,5", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("4\\/5\\/6", "3,5,6", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("1\\/2\\/6", "1,2,6", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("\\/5\\/6", "5,6", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("\\/4\\/5", "4,5", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("5\\/6", "5,6", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("3\\/5", "3,5", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("2\\/3", "2,3", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("1\\/2", "1,2", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("4\\/6", "4,6", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("3\\/4", "3,4", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("/2", "2", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("/4", "4", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("/6", "6", layer_marker_Zeng)) %>%
  mutate(layer_marker_Zeng = gsub("\\/", NA, layer_marker_Zeng))

layer_marker_per_source <- full_join(He_layer_marker, Maynard_layer_enrichment, 
                                     by = 'gene_symbol',
                                     suffix = c('_He', '_Maynard')) %>%
  full_join(Zeng_marker_annotation, by = "gene_symbol") %>%
  mutate(layer_marker_He = gsub("Layer", "", layer_marker_He)) %>%
  mutate(layer_marker_Maynard = gsub("Layer", "", layer_marker_Maynard))

write.csv(layer_marker_per_source, here("Data", "processed_data",
                                        ""))
