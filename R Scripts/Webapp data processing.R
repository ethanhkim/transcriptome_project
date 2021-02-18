## Add WM NaN values to MTG data
AIBS_MTG_df <- read.csv(here("Data", "Allen", "MTG_df_01_21.csv")) %>%
  dplyr::select(-X, -X1)

# Filter by cell type
filter_by_cell_type <- function(AIBS_df, cell_type) {
  df <- AIBS_df
  
  df %<>%
    filter(class_label == cell_type)
  
  return(df)
}

# Z-score and add NaN to each filtered df
add_nan <- function(filtered_df) {
  df <- filtered_df
  
  df %<>% pivot_longer(names_from = cortical_layer_label,
                       values_from = mean_expression) %>%
    add_column
}

test <- AIBS_MTG_df %>%
  filter(class_label == "GABAergic")

AIBS_MTG_GABA <- AIBS_MTG_df %>%
  filter_by_cell_type(cell_type = "GABAergic") %>%
  pivot_wider(names_from = cortical_layer_label,
              values_from = mean_expression)

test <- AIBS_MTG_df %>%
  filter_by_cell_type(cell_type = "GABAergic") %>%
  spread(cortical_layer_label, mean_expression)
AIBS_MTG_GLUT <- AIBS_MTG_df %>%
  filter_by_cell_type(cell_type = "Glutamatergic")
AIBS_MTG_NONN <- AIBS_MTG_df %>%
  filter_by_cell_type(cell_type = "Non-neuronal")

AIBS_MTG_df %<>% pivot_wider(names_from = cortical_layer_label,
                             values_from = mean_expression)
  dplyr::select(-cortical_layer_label) %>%
  add_column(WM = NaN) %>%
  pivot_longer(cols = L1:WM, names_to = 'test', values_to = 'mean-expression')
  
