#Create merged data between human datasets: Zeng et al, He et al, Maynard et al

#Create common gene list between Zeng, He, Maynard

commongene_List_he_maynard_zeng <- intersect(sce_layer_data$gene_symbol, common_geneList) %>%
  na.omit() %>%
  as_vector()

#select out rows containing only gene symbol from above lists
filtered_He <- filtered_He[filtered_He$gene_symbol %in% commongene_List_he_maynard_zeng,]
filtered_Zeng <- updatedZengTable[updatedZengTable$gene_symbol %in% commongene_List_he_maynard_zeng,]
filtered_Maynard <- sce_layer_data[sce_layer_data$gene_symbol %in% commongene_List_he_maynard_zeng,] 
filtered_Maynard <- filtered_Maynard[-c(265, 785),]


Zeng_Layer1 <- filter(filtered_Zeng_2, Cortical.marker..human. == "layer 1")
Zeng_Layer2 <- filter(filtered_Zeng_2, Cortical.marker..human. == "layer 2")
Zeng_Layer3 <- filter(filtered_Zeng_2, Cortical.marker..human. == "layer 3")
Zeng_Layer4 <- filter(filtered_Zeng_2, Cortical.marker..human. == "layer 4")
Zeng_Layer5 <- filter(filtered_Zeng_2, Cortical.marker..human. == "layer 5")
Zeng_Layer6 <- filter(filtered_Zeng_2, Cortical.marker..human. == "layer 6")

#Merge layer data by gene symbol

merged_layer1 <- merge(He_Layer1, Zeng_Layer1, by = "gene_symbol", all.x = TRUE)
merged_layer2 <- merge(He_Layer2, Zeng_Layer2, by = "gene_symbol", all.x = TRUE)
merged_layer3 <- merge(He_Layer3, Zeng_Layer3, by = "gene_symbol", all.x = TRUE)
merged_layer4 <- merge(He_Layer4, Zeng_Layer4, by = "gene_symbol", all.x = TRUE)
merged_layer5 <- merge(He_Layer5, Zeng_Layer5, by = "gene_symbol", all.x = TRUE)
merged_layer6 <- merge(He_Layer6, Zeng_Layer6, by = "gene_symbol", all.x = TRUE)

