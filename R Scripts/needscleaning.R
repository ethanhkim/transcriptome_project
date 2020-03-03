#scripts that need work/cleaning to optimize

if (!("V1") %in% column_names) {
  layer1 <- filtered_DS1_Human %>%
    dplyr::select(gene_symbol, V1)
} else if (("V2") %in% column_names) {
  layer2 <- filtered_DS1_Human %>%
    dplyr::select(gene_symbol, V2, V3, V4)
} else if (("V4" | "V5" | "V6") %in% column_names) {
  layer3 <- filtered_DS1_Human %>%
    dplyr::select(gene_symbol, V4, V5, V6) 
} else if (("V7" | "V8" | "V9") %in% column_names) {
  layer4 <- filtered_DS1_Human %>%
    dplyr::select(gene_symbol, V7, V8, V9)
} 