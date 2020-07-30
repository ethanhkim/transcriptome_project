#Create merged data between human datasets: Zeng et al, He et al, Maynard et al

#Create common gene list between =He, Maynard

common_geneList_he_zeng <- intersect(maynard_geneList, he_geneList)

#select out rows containing only gene symbol from above lists
filtered_He <- DS1_Human[DS1_Human$gene_symbol %in% common_geneList_he_zeng,]
rownames(filtered_He) <- c(1:18271)
filtered_Maynard <- sce_layer_data[sce_layer_data$gene_symbol %in% common_geneList_he_zeng,] 
filtered_Maynard$row <- c(1:18266)
filtered_Maynard <- filtered_Maynard[, c(78,79,1:77)]

#Remove duplicates#
#He et al duplicates#
unique(filtered_He$gene_symbol[duplicated(filtered_He$gene_symbol) | duplicated(filtered_He$gene_symbol, fromLast=TRUE)])

filtered_He <- filtered_He[-c(16526, 17875, 2663, 18215, 3227, 13404, 16693, 18078, 16412, 5855, 18154, 12445, 16757, 
                              17199, 17051, 17780, 17713, 17999, 17166, 17453, 17803, 17233, 11376, 15827, 15826, 15768, 
                              16550, 16547, 18046, 17017, 18106, 15626, 17703, 18210, 16968, 18182, 17802, 17865, 17843, 9441),]

#Maynard et al duplicates#
unique(filtered_Maynard$gene_symbol[duplicated(filtered_Maynard$gene_symbol) | duplicated(filtered_Maynard$gene_symbol, fromLast=TRUE)])

filtered_Maynard <- filtered_Maynard[-c(5442, 1352, 2619, 2675, 3226, 3446, 3540, 3937, 5059, 5060, 5346, 54142, 5959, 6416, 
                                        6504, 6517, 6772, 7149, 7847, 8094, 8762, 9048, 9118, 9242, 10576, 10816, 11360, 
                                        12116, 12236, 13354, 13739, 14024, 15121, 15126, 16894, 18083),]
filtered_Maynard <- filtered_Maynard[,-c(2, 79)]

## Separate by layer ##

# He
He_Layer1 <- filtered_He %>% dplyr::select(gene_symbol, V1)
He_Layer2 <- filtered_He %>% dplyr::select(gene_symbol, V2, V3)
He_Layer3 <- filtered_He %>% dplyr::select(gene_symbol, V5, V6)
He_Layer4 <- filtered_He %>% dplyr::select(gene_symbol, V7, V8)
He_Layer5 <- filtered_He %>% dplyr::select(gene_symbol, V10, V11, V12)
He_Layer6 <- filtered_He %>% dplyr::select(gene_symbol, V13, V14, V15, V16)

# Maynard

Maynard_151507 <- filtered_Maynard[,c(1, 2:8)]
Maynard_151508 <- filtered_Maynard[,c(1, 9:15)]
Maynard_151509 <- filtered_Maynard[,c(1, 16:22)]
Maynard_151510 <- filtered_Maynard[,c(1, 23:29)]
Maynard_151669 <- filtered_Maynard[,c(1, 30:34)]
Maynard_151670 <- filtered_Maynard[,c(1, 35:39)]
Maynard_151671 <- filtered_Maynard[,c(1, 40:44)]
Maynard_151672 <- filtered_Maynard[,c(1, 45:49)]
Maynard_151673 <- filtered_Maynard[,c(1, 50:56)]
Maynard_151674 <- filtered_Maynard[,c(1, 57:63)]
Maynard_151675 <- filtered_Maynard[,c(1, 64:70)]
Maynard_151676 <- filtered_Maynard[,c(1, 71:77)]
