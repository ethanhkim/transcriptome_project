## AIBS normalized data ##

# Load necessary libraries #

library(data.table)
library(org.Hs.eg.db)
library(dplyr)
library(magrittr)
library(stringr)
library(conflicted)
library(biomaRt)
library(HGNChelper)
library(edgeR)
library(here)

# Set conflicts #
conflict_prefer('select', 'dplyr')
conflict_prefer('filter', 'dplyr')
conflict_prefer('cpm', 'edgeR')
conflict_prefer('rename', 'dplyr')


# Function to separate out by cell type

separate_by_cell_type <- function(AIBS_MTG_df, cell_type) {
  df <- AIBS_MTG_df

  df %<>% 
    # Filter by specified cell type
    filter(class_label == cell_type) %>%
    # Widen to gene by layer
    pivot_wider(names_from = cortical_layer_label,
                values_from = mean_expression) %>%
    # Gene symbols to row name
    column_to_rownames(var = 'gene_symbol')
  
  return(df)
}


AIBS_MTG <- fread(here('Data', 'Allen', 'MTG_df_01_21.csv')) %>%
  # remove col number column
  select(-V1)

AIBS_MTG_cell_type <- list()
for (i in c("GABAergic", "Glutamatergic", "Non-neuronal")) {
  AIBS_MTG_cell_type[[i]] <- separate_by_cell_type(AIBS_MTG, i)
}

AIBS_MTG_all_cell_types <- rbind(AIBS_MTG_cell_type$GABAergic, 
                               AIBS_MTG_cell_type$Glutamatergic, 
                               AIBS_MTG_cell_type$`Non-neuronal`)




