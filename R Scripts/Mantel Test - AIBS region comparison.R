## Mantel tests for between scRNA-seq data ##

# Load AIBS_logCPM data
load(here("Data", "processed_data", "AIBS_logCPM.RData"))

# Separate by region and convert long to wide

# Primary auditory
A1C_df <- AIBS_logCPM %>%
  filter(region_label == "A1C") %>%
  group_by(gene_symbol, cortical_layer_label) %>%
  summarize(mean_exp = mean(log_mean_expression)) %>%
  spread(cortical_layer_label, mean_exp) %>%
  mutate(L6 = (L6a+L6b)/2) %>%
  
  
# Convert long to wide format
CgG_df %<>%
  group_by(gene_symbol, cortical_layer_label) %>%
  summarize(mean_exp = mean(log_mean_expression)) %>%
  spread(cortical_layer_label, mean_exp) %>%
  mutate(L5 = (L5a+L5b)/2) %>%
  add_column(L3 = 0, L4 = 0, WM = 0) %>%
  select(gene_symbol, L1, L2, L3, L4, L5, L6)