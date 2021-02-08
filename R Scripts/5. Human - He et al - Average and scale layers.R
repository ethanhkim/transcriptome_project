## Load required libraries ##

library(tidyverse)
library(magrittr)
library(conflicted)

conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")

## Average He et al ##

He_Layer1_averaged <- He_DS1_Human %>% 
  dplyr::select(gene_symbol, V1)

He_Layer2_averaged <- He_DS1_Human %>% 
  dplyr::select(gene_symbol, V2, V3, V4) %>%
  mutate(V4_weighted = V4 * 0.5) %>%
  rowwise() %>%
  mutate(wt.mean = mean(c(V2, V3, V4_weighted)))

He_Layer3_averaged <- He_DS1_Human %>% 
  select(gene_symbol, V4, V5, V6) %>%
  mutate(V4_weighted = V4 * 0.5) %>%
  rowwise() %>%
  mutate(wt.mean = mean(c(V4_weighted, V5, V6)))

He_Layer4_averaged <- He_DS1_Human %>% 
  select(gene_symbol, V7, V8, V9) %>%
  mutate(V9_weighted = V9 * 0.5) %>%
  rowwise() %>%
  mutate(wt.mean = mean(c(V7, V8, V9_weighted)))

He_Layer5_averaged <- He_DS1_Human %>% 
  select(gene_symbol, V9, V10, V11, V12) %>%
  mutate(V9_weighted = V9 * 0.5) %>%
  rowwise() %>%
  mutate(wt.mean = mean(c(V9_weighted, V10, V11, V12)))

He_Layer6_averaged <- He_DS1_Human %>% 
  select(gene_symbol, V13, V14, V15, V16) %>%
  rowwise() %>%
  mutate(mean = mean(c(V13, V14, V15, V16)))

He_WM_averaged <- He_DS1_Human  %>%
  select(gene_symbol, V17)


## Create averaged He et al. data table from the scaled values
He_DS1_averaged_by_layer <- tibble(
  gene_symbol = He_DS1_Human$gene_symbol,
  Layer_1 = He_Layer1_averaged$V1,
  Layer_2 = He_Layer2_averaged$wt.mean,
  Layer_3 = He_Layer3_averaged$wt.mean,
  Layer_4 = He_Layer4_averaged$wt.mean,
  Layer_5 = He_Layer5_averaged$wt.mean,
  Layer_6 = He_Layer6_averaged$mean,
  Layer_7 = He_WM_averaged$V17
)


## Cleanup workspace
rm(He_Layer1_averaged, He_Layer2_averaged, He_Layer3_averaged, He_Layer4_averaged, He_Layer5_averaged, He_Layer6_averaged, He_WM_averaged)

## Scale Layers ##

He_values_scaled <- He_values %>%
  dplyr::select(-"gene_symbol") %>%
  t() %>%
  scale() %>%
  t() %>%
  as_tibble() %>%
  add_column(He_values$gene_symbol) %>%
  rename(gene_symbol= "He_values$gene_symbol") %>%
  dplyr::select("gene_symbol", everything())


## Transpose data ##

#Transpose table such that table represents the 17 separate cuts for each gene
rownames(He_values) <- common_gene_list 
He_values_transposed <- He_values %>%
  dplyr::select(-"gene_symbol") %>%
  t() %>%
  as_tibble()



