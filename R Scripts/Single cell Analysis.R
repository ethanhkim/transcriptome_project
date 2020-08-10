library(tidyverse)
library(here)

allen_singlecellMetadata <- read_csv(here("Data", "Allen", "metadata.csv"))
View(allen_singlecellMetadata)

unique(allen_singlecellMetadata$cortical_layer_color)
