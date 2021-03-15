## Load required libraries ##

library(tidyverse)
library(magrittr)
library(conflicted)

conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")

## Average He et al ##

He_Layer1_averaged <- He_DS1_Human %>% 
  dplyr::select(gene_symbol, S1)

He_Layer2_averaged <- He_DS1_Human %>% 
  dplyr::select(gene_symbol, S2, S3, S4) %>%
  mutate(S4_weighted = S4 * 0.5) %>%
  rowwise() %>%
  mutate(wt.mean = mean(c(S2, S3, S4_weighted)))

He_Layer3_averaged <- He_DS1_Human %>% 
  select(gene_symbol, S4, S5, S6) %>%
  mutate(S4_weighted = S4 * 0.5) %>%
  rowwise() %>%
  mutate(wt.mean = mean(c(S4_weighted, S5, S6)))

He_Layer4_averaged <- He_DS1_Human %>% 
  select(gene_symbol, S7, S8, S9) %>%
  mutate(S9_weighted = S9 * 0.5) %>%
  rowwise() %>%
  mutate(wt.mean = mean(c(S7, S8, S9_weighted)))

He_Layer5_averaged <- He_DS1_Human %>% 
  select(gene_symbol, S9, S10, S11, S12) %>%
  mutate(S9_weighted = S9 * 0.5) %>%
  rowwise() %>%
  mutate(wt.mean = mean(c(S9_weighted, S10, S11, S12)))

He_Layer6_averaged <- He_DS1_Human %>% 
  select(gene_symbol, S13, S14, S15, S16) %>%
  rowwise() %>%
  mutate(mean = mean(c(S13, S14, S15, S16)))

He_WM_averaged <- He_DS1_Human  %>%
  select(gene_symbol, S17)


## Create averaged He et al. data table from the scaled values
He_DS1_averaged_by_layer <- tibble(
  gene_symbol = He_DS1_Human$gene_symbol,
  Layer_1 = He_Layer1_averaged$S1,
  Layer_2 = He_Layer2_averaged$wt.mean,
  Layer_3 = He_Layer3_averaged$wt.mean,
  Layer_4 = He_Layer4_averaged$wt.mean,
  Layer_5 = He_Layer5_averaged$wt.mean,
  Layer_6 = He_Layer6_averaged$mean,
  WM = He_WM_averaged$S17
)


## Cleanup workspace
rm(He_Layer1_averaged, He_Layer2_averaged, He_Layer3_averaged, He_Layer4_averaged, 
   He_Layer5_averaged, He_Layer6_averaged, He_WM_averaged, He_DS1_Human)

save(He_DS1_averaged_by_layer, 
     file = here("Data", "processed_data", "He_DS1_averaged_by_layer.Rdata"))

### OUTDATED ###

## Scale Layers ##

#He_values_scaled <- He_values %>%
#  dplyr::select(-"gene_symbol") %>%
#  t() %>%
#  scale() %>%
#  t() %>%
#  as_tibble() %>%
#  add_column(He_values$gene_symbol) %>%
#  rename(gene_symbol= "He_values$gene_symbol") %>%
#  dplyr::select("gene_symbol", everything())


## Transpose data ##

#Transpose table such that table represents the 17 separate cuts for each gene
#rownames(He_values) <- common_gene_list 
#He_values_transposed <- He_values %>%
#  dplyr::select(-"gene_symbol") %>%
#  t() %>%
#  as_tibble()

