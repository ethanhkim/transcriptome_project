#Unify gene_symbol column ID

cleanedZengTable <- dplyr::rename(cleanedZengTable, gene_symbol = "Gene.symbol")

#Create gene lists containing unique values

heList <- DS1_Human$gene_symbol
zengList <- unique(cleanedZengTable$gene_symbol)

#Create list of common genes
commonList <- dplyr::intersect(heList, zengList)

#select out rows containing only genes from above lists
filtered_He <- DS1_Human[DS1_Human$gene_symbol %in% commonList,]
filtered_Zeng <- cleanedZengTable[cleanedZengTable$gene_symbol %in% commonList,]

#Note - SERF1A, TUT1, CRHR1 contain duplicates

#Separate by layer 

He_Layer1 <- filtered_He %>% dplyr::select(gene_symbol, V1)
He_Layer2 <- filtered_He %>% dplyr::select(gene_symbol, V2, V3, V4)
He_Layer3 <- filtered_He %>% dplyr::select(gene_symbol, V4, V5, V6)
He_Layer4 <- filtered_He %>% dplyr::select(gene_symbol, V7, V8, V9)
He_Layer5 <- filtered_He %>% dplyr::select(gene_symbol, V9, V10, V11, V12)
He_Layer6 <- filtered_He %>% dplyr::select(gene_symbol, V13, V14, V15, V16)

Zeng_Layer1 <- filter(filtered_Zeng, Cortical.marker..human. == "layer 1")
Zeng_Layer2 <- filter(filtered_Zeng, Cortical.marker..human. == "layer 2")
Zeng_Layer3 <- filter(filtered_Zeng, Cortical.marker..human. == "layer 3")
Zeng_Layer4 <- filter(filtered_Zeng, Cortical.marker..human. == "layer 4")
Zeng_Layer5 <- filter(filtered_Zeng, Cortical.marker..human. == "layer 5")
Zeng_Layer6 <- filter(filtered_Zeng, Cortical.marker..human. == "layer 6")

#Merge layer data by gene symbol

merged_layer1 <- merge(He_Layer1, Zeng_Layer1, by = "gene_symbol", all.x = TRUE)
merged_layer2 <- merge(He_Layer2, Zeng_Layer2, by = "gene_symbol", all.x = TRUE)
merged_layer3 <- merge(He_Layer3, Zeng_Layer3, by = "gene_symbol", all.x = TRUE)
merged_layer4 <- merge(He_Layer4, Zeng_Layer4, by = "gene_symbol", all.x = TRUE)
merged_layer5 <- merge(He_Layer5, Zeng_Layer5, by = "gene_symbol", all.x = TRUE)
merged_layer6 <- merge(He_Layer6, Zeng_Layer6, by = "gene_symbol", all.x = TRUE)
