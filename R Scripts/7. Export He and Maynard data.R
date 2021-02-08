## Export data to Webapp ##

#Maynard layer enrichment code
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

Maynard_dataset_average <- Maynard_dataset_average %>%
  arrange(gene_symbol) %>%
  filter(gene_symbol %in% Maynard_layer_enrichment$gene_symbol)

Maynard_layer_enrichment %<>% filter(gene_symbol %in% Maynard_dataset_average$gene_symbol)

Maynard_dataset_average %<>%
  column_to_rownames(var = "gene_symbol") %>%
  t() %>%
  scale() %>%
  t() %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "gene_symbol") %>%
  add_column(
    marker_label = Maynard_layer_enrichment$layer_marker_label,
    Layer_1_marker = Maynard_layer_enrichment$layer_marker_label_1,
    Layer_2_marker = Maynard_layer_enrichment$layer_marker_label_2,
    Layer_3_marker = Maynard_layer_enrichment$layer_marker_label_3,
    Layer_4_marker = Maynard_layer_enrichment$layer_marker_label_4,
    Layer_5_marker = Maynard_layer_enrichment$layer_marker_label_5,
    Layer_6_marker = Maynard_layer_enrichment$layer_marker_label_6,
    Layer_WM_marker = Maynard_layer_enrichment$layer_marker_label_WM)

He_DS1_averaged_by_layer_export <- He_DS1_averaged_by_layer %>%
  column_to_rownames(var = "gene_symbol") %>%
  t() %>%
  scale() %>%
  t() %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "gene_symbol")

## He Layer markers ##
He_LM_Path <- here("Data", "He et al", "Supplementary Table 2.xlsx")
He_Layer_markers <- read_xlsx(path = He_LM_Path) %>%
  dplyr::select("Gene symbol", "Layer marker in human") %>%
  mutate_at(vars("Layer marker in human"), na_if, "No") %>%
  mutate_at(vars("Layer marker in human"), na_if, "NA") %>%
  dplyr::rename(gene_symbol = "Gene symbol",  layer_marker = "Layer marker in human")

He_DS1_Human_averaged <-  merge(x = He_DS1_averaged_by_layer_export, y = He_Layer_markers, 
                                       by = "gene_symbol", all.x = TRUE) %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  rename(Layer_WM = "Layer_7") %>%
  mutate(layer_marker_label = ifelse(test = (is.na(layer_marker)),
                                     yes = NA,
                                     no = "*")) %>%
  mutate(Layer_1_marker = ifelse(test = (layer_marker == "L1"),
                                       yes = "Layer_1",
                                       no = NA)) %>%
  mutate(Layer_2_marker = ifelse(test = (layer_marker == "L2"),
                                       yes = "Layer_2",
                                       no = NA)) %>%
  mutate(Layer_3_marker = ifelse(test = (layer_marker == "L3"),
                                       yes = "Layer_3",
                                       no = NA)) %>%
  mutate(Layer_4_marker = ifelse(test = (layer_marker == "L4"),
                                       yes = "Layer_4",
                                       no = NA)) %>%
  mutate(Layer_5_marker = ifelse(test = (layer_marker == "L5"),
                                       yes = "Layer_5",
                                       no = NA)) %>%
  mutate(Layer_6_marker = ifelse(test = (layer_marker == "L6"),
                                       yes = "Layer_6",
                                       no = NA)) %>%
  mutate(Layer_WM_marker = ifelse(test = (layer_marker == "WM"),
                                   yes = "WM",
                                   no = NA)) %>%
  rename(marker_label = "layer_marker_label", WM = "Layer_WM") %>%
  dplyr::select(-"layer_marker")



write_csv(He_DS1_Human_averaged, './R Scripts/export_data/He_DS1_Human_averaged.csv')
save(He_DS1_Human_averaged, file = './R Scripts/export_data/He_DS1_Human_averaged.Rdata')


write_csv(Maynard_dataset_average, './R Scripts/export_data/Maynard_dataset_average.csv')
save(Maynard_dataset_average, file = './R Scripts/export_data/Maynard_dataset_average.Rdata')

                                       
  