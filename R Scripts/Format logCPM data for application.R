## Format He, Maynard ad Allen data for webapp

# Load libraries
library(dplyr)
library(tidyr)
library(tibble)
library(here)

# Load data
load(here("Data", "processed_data", "MTG_logCPM_dataset.Rdata"))
load(here("Data", "processed_data", "He_DS1_logCPM_dataset.Rdata"))
load(here("Data", "processed_data", "Maynard_logCPM_dataset.Rdata"))

# Z-score log-transformed data

# Function to add NaN values to "white matter"
add_NaN_to_AIBS <- function(AIBS_MTG_df, cell_type) {
  
  df <- AIBS_MTG_df %>%
    filter(class_label == cell_type) %>%
    add_column(WM = NaN) %>%
    pivot_longer(cols = L1:WM,
                 names_to = "cortical_layer_label",
                 values_to = "expression") 
  
  return(df)
}

# Create a list of dataframes of lengthened
# Allen data with NaN white matter values inserted
AIBS_MTG_cell_type <- list()
for (i in c("GABAergic", "Glutamatergic", "Non-neuronal")) {
  AIBS_MTG_cell_type[[i]] <- add_NaN_to_AIBS(MTG_logCPM_dataset, i)
}

MTG_logCPM_dataset_lengthened <- rbind(AIBS_MTG_cell_type$GABAergic, 
                                       AIBS_MTG_cell_type$Glutamatergic, 
                                       AIBS_MTG_cell_type$`Non-neuronal`)

save(MTG_logCPM_dataset_lengthened, 
     file = here("Data", "processed_data", 
                 "MTG_logCPM_dataset_lengthened.Rdata"))

## He et al data



test_AIBS <- MTG_logCPM_dataset %>%
  pivot_longer(cols = L1:L6,
               names_to = "cortical_layer_label",
               values_to = "expression")

min(test_AIBS$expression) #-11.22733
max(test_AIBS$expression) #15.50301
hist(test_AIBS$expression)

test_Maynard <- Maynard_logCPM_dataset %>%
  pivot_longer(cols = Layer_1:WM,
               names_to = "cortical_layer_label",
               values_to = "expression")

min(test_Maynard$expression) #-2.671885
max(test_Maynard$expression) #14.14513
hist(test_Maynard$expression)


test_He <- He_DS1_logCPM_dataset %>%
  pivot_longer(cols = Layer_1:WM,
               names_to = "cortical_layer_label",
               values_to = "expression")


min(test_He$expression) #-5.41005
max(test_He$expression) #15.56536
hist(test_He$expression)


  
