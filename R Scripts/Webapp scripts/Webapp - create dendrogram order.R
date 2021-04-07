library(ggplot2)
library(magrittr)
library(plotly)
library(tibble)
library(shinyjs)
library(tidyr)
library(dplyr)
library(purrr)
library(stringr)
library(stringi)
library(here)
library(scales)
library(ggdendro)
library(data.table)

conflict_prefer("filter", "dplyr")

ordered_data <- He_DS1_logCPM_dataset %>%
  replace_na(list(
    Layer_1 = 0,
    Layer_2 = 0,
    Layer_3 = 0,
    Layer_4 = 0,
    Layer_5 = 0,
    Layer_6 = 0,
    WM = 0
  )) %>%
  mutate_if(is.numeric,as.character, is.factor, as.character) %>%
  select(gene_symbol:WM)

ordered_data_matrix <- as.matrix(test_ordered_data)
rownames(ordered_data) <- ordered_data$gene_symbol
data_dendro <- as.dendrogram(hclust(d = dist(x = ordered_data_matrix)))

data_order <- order.dendrogram(data_dendro)

save(data_order, file = here("Data", "processed_data", "He_ordered_index.Rdata"))

