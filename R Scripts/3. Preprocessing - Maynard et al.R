## Maynard Data ##

library(spatialLIBD)
library(tidyverse)
library(org.Hs.eg.db)
library(here)

sce_layer <- fetch_data(type = "sce_layer")
Maynard_dataset <- as_tibble(sce_layer@assays@data@listData$logcounts)
Maynard_ensembl_list <- sce_layer@rowRanges@ranges@NAMES
Maynard_dataset$Ensembl_ID <- Maynard_ensembl_list



#Change Ensembl ID to gene symbol
Maynard_dataset$gene_symbol <- mapIds(org.Hs.eg.db, keys = Maynard_dataset$Ensembl_ID, keytype = "ENSEMBL", column="SYMBOL")
Maynard_dataset <- Maynard_dataset %>%
  dplyr::select(-"Ensembl_ID") %>% 
  dplyr::select("gene_symbol", everything()) %>%
  filter(!is.na(gene_symbol)) %>%
  distinct(gene_symbol, .keep_all = TRUE)

# Write raw layer-level dataset for Maynard
# unable to install spatialLIBD on SCC.

write.csv(Maynard_dataset, here('Data', 'raw_data', 'Maynard et al',
                                'layer_level_data.csv'),
          row.names = FALSE)

## Separate by layers and add average column
Maynard_dataset_Layer1 <- Maynard_dataset %>%
  dplyr::select("gene_symbol", contains("Layer1")) %>%
  mutate(dataset_mean = rowMeans(dplyr::select(., `151507_Layer1`:`151676_Layer1`)))

Maynard_dataset_Layer2 <- Maynard_dataset %>%
  dplyr::select("gene_symbol", contains("Layer2")) %>%
  mutate(dataset_mean = rowMeans(dplyr::select(., `151507_Layer2`:`151676_Layer2`)))

Maynard_dataset_Layer3 <- Maynard_dataset %>%
  dplyr::select("gene_symbol", contains("Layer3")) %>%
  dplyr::select(-"151669_Layer3", -"151670_Layer3", -"151671_Layer3", -"151672_Layer3") %>%
  mutate(dataset_mean = rowMeans(dplyr::select(., `151507_Layer3`:`151676_Layer3`)))

Maynard_dataset_Layer4 <- Maynard_dataset %>%
  dplyr::select("gene_symbol", contains("Layer4")) %>%
  dplyr::select(-"151669_Layer4", -"151670_Layer4", -"151671_Layer4", -"151672_Layer4") %>%
  mutate(dataset_mean = rowMeans(dplyr::select(., `151507_Layer4`:`151676_Layer4`)))

Maynard_dataset_Layer5 <- Maynard_dataset %>%
  dplyr::select("gene_symbol", contains("Layer5")) %>%
  dplyr::select(-"151669_Layer5", -"151670_Layer5", -"151671_Layer5", -"151672_Layer5") %>%
  mutate(dataset_mean = rowMeans(dplyr::select(., `151507_Layer5`:`151676_Layer5`)))

Maynard_dataset_Layer6 <- Maynard_dataset %>%
  dplyr::select("gene_symbol", contains("Layer6")) %>%
  dplyr::select(-"151669_Layer6", -"151670_Layer6", -"151671_Layer6", -"151672_Layer6") %>%
  mutate(dataset_mean = rowMeans(dplyr::select(., `151507_Layer6`:`151676_Layer6`)))

Maynard_dataset_WM <- Maynard_dataset %>%
  dplyr::select("gene_symbol", contains("WM")) %>%
  dplyr::select(-"151669_WM", -"151670_WM", -"151671_WM", -"151672_WM") %>%
  mutate(dataset_mean = rowMeans(dplyr::select(., `151507_WM`:`151676_WM`)))


Maynard_dataset_average <- tibble(
  gene_symbol = Maynard_dataset$gene_symbol,
  Layer_1 = Maynard_dataset_Layer1$dataset_mean,
  Layer_2 = Maynard_dataset_Layer2$dataset_mean,
  Layer_3 = Maynard_dataset_Layer3$dataset_mean,
  Layer_4 = Maynard_dataset_Layer4$dataset_mean,
  Layer_5 = Maynard_dataset_Layer5$dataset_mean,
  Layer_6 = Maynard_dataset_Layer6$dataset_mean,
  WM = Maynard_dataset_WM$dataset_mean,
)

# Clean up workspace
rm(Maynard_dataset_Layer1, Maynard_dataset_Layer2, Maynard_dataset_Layer3, Maynard_dataset_Layer4, Maynard_dataset_Layer5,
   Maynard_dataset_Layer6, Maynard_dataset_WM)

## Separate by sample ID
#Maynard_151507 <- Maynard_dataset %>%
#  dplyr::select("gene_symbol", contains("151507"))
#Maynard_151508 <- Maynard_dataset %>%
#  dplyr::select("gene_symbol", contains("151508"))
#Maynard_151509 <- Maynard_dataset %>%
#  dplyr::select("gene_symbol", contains("151509"))
#Maynard_151510 <- Maynard_dataset %>%
#  dplyr::select("gene_symbol", contains("151510"))
#Maynard_151669 <- Maynard_dataset %>%
#  dplyr::select("gene_symbol", contains("151669"))
#Maynard_151670 <- Maynard_dataset %>%
#  dplyr::select("gene_symbol", contains("151670"))
#Maynard_151671 <- Maynard_dataset %>%
#  dplyr::select("gene_symbol", contains("151671"))
#Maynard_151672 <- Maynard_dataset %>%
#  dplyr::select("gene_symbol", contains("151672"))
#Maynard_151673 <- Maynard_dataset %>%
#  dplyr::select("gene_symbol", contains("151673"))
#Maynard_151674 <- Maynard_dataset %>%
#  dplyr::select("gene_symbol", contains("151674"))
#Maynard_151675 <- Maynard_dataset %>%
#  dplyr::select("gene_symbol", contains("151675"))
#Maynard_151676 <- Maynard_dataset %>%
#  dplyr::select("gene_symbol", contains("151676"))


