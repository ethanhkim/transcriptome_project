library(tidyverse)
library(reshape2)

He_DS1_filtered <- He_DS1_averaged_by_layer %>%
  dplyr::slice(match(He_Maynard_common_gene_list, gene_symbol)) %>%
  column_to_rownames(var = "gene_symbol") %>%
  t()


Maynard_filtered <- Maynard_dataset_average %>%
  dplyr::slice(match(He_Maynard_common_gene_list, gene_symbol)) %>%
  column_to_rownames(var = "gene_symbol") %>%
  t()

test <- cor(He_DS1_filtered, Maynard_151507_filtered, method = "pearson") %>%
  melt() %>%
  as_tibble()

test %>%
  group_by(Var1 == Var2) %>%
  summarise(mean_corr = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE)) 





