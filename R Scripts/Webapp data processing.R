## Add WM NaN values to MTG data

AIBS_MTG_df <- fread(here("Data", "Allen", "MTG_df_01_21.csv")) %>%
  select(-V1, -X1)

add_NaN_to_AIBS <- function(AIBS_MTG_df, cell_type) {
  
  df <- AIBS_MTG_df %>%
    filter(class_label == cell_type) %>%
    pivot_wider(names_from = cortical_layer_label,
                values_from = mean_expression,
                values_fn = list)
  
  listcol_to_col <- function(column) {
    mutated_col <- do.call(rbind.data.frame, column) %>%
      pull(var = 1)
    return(mutated_col)
  }
  
  mutated_df <- data.frame(df[1:2], 
                   sapply(df[3:8], listcol_to_col)) %>%
    add_column(WM = NaN) %>%
    pivot_longer(cols = L1:WM,
                 names_to = "cortical_layer_label",
                 values_to = "mean_expression")
  
  return(mutated_df)
}


AIBS_MTG_cell_type <- list()
for (i in c("GABAergic", "Glutamatergic", "Non-neuronal")) {
  AIBS_MTG_cell_type[[i]] <- add_NaN_to_AIBS(AIBS_MTG_df, i)
}

AIBS_logCPM_dataset <- rbind(AIBS_MTG_cell_type$GABAergic, 
                      AIBS_MTG_cell_type$Glutamatergic, 
                      AIBS_MTG_cell_type$`Non-neuronal`)

save(AIBS_logCPM_dataset, file = here("Data", "processed_data", "AIBS_logCPM_dataset.Rdata"))


  
