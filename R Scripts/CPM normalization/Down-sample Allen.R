GABA <- MTG %>%
  filter(class_label == 'GABAergic')

GLUT <- MTG %>%
  filter(class_label == 'Glutamatergic')

NONN <- MTG %>%
  filter(class_label == 'Non-neuronal')

GABA_layer_nuc_count <- GABA %>%
  group_by(cortical_layer_label) %>%
  count()

GLUT_layer_nuc_count <- GLUT %>%
  group_by(cortical_layer_label) %>%
  count()

NONN_layer_nuc_count <- NONN %>%
  group_by(cortical_layer_label) %>%
  count()

cell_type_nuc_count <- tibble(
  GABAergic = GABA_layer_nuc_count %>% pull(),
  Glutamatergic = GLUT_layer_nuc_count %>% pull(),
  Non_neuronal = NONN_layer_nuc_count %>% pull()
)

sample_nuclei <- function(df, layer_label, n_sample) {
  
  df %<>% filter(cortical_layer_label == layer_label)
  
  sampled_layer <- df[sample(nrow(df), n_sample),]
  
  return(sampled_layer)
}


GABA_sampled <- list()
GLUT_sampled <- list()
NONN_sampled <- list()
for (layer in c("L1", "L2", "L3", "L4", "L5", "L6")) {
  GABA_sampled[[layer]] <- sample_nuclei(GABA, layer, 381)
  GLUT_sampled[[layer]] <- sample_nuclei(GLUT, layer, 274)
  NONN_sampled[[layer]] <- sample_nuclei(NONN, layer, 125)
}

GABA_sampled_df <- rbindlist(GABA_sampled)
GLUT_sampled_df <- rbindlist(GLUT_sampled)
NONN_sampled_df <- rbindlist(NONN_sampled)

AIBS_sampled_df <- rbind(GABA_sampled_df, GLUT_sampled_df, NONN_sampled_df)


AIBS_sampled_cell_sum_count <- sum_gene_count(AIBS_sampled_df, "gene")

write.csv(AIBS_sampled_cell_sum_count, 
          file = here("Data", "raw_data", "Allen", "AIBS_sampled_cell_sum_count.csv"))
