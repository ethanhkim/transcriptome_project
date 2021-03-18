# Zeng layer marker table

Zeng_L1 <- read.delim(here("Data", "genelists", "Zeng.Layer.1.Markers.txt"), header=F) %>%
  as_tibble() %>%
  add_column(layer = 1)
Zeng_L2 <- read.delim(here("Data", "genelists", "Zeng.Layer.2.Markers.txt"), header=F) %>%
  as_tibble() %>%
  add_column(layer = 2)
Zeng_L3 <- read.delim(here("Data", "genelists", "Zeng.Layer.3.Markers.txt"), header=F) %>%
  as_tibble() %>%
  add_column(layer = 3)
Zeng_L4 <- read.delim(here("Data", "genelists", "Zeng.Layer.4.Markers.txt"), header=F) %>%
  as_tibble() %>%
  add_column(layer = 4)
Zeng_L5 <- read.delim(here("Data", "genelists", "Zeng.Layer.5.Markers.txt"), header=F) %>%
  as_tibble() %>%
  add_column(layer = 5)
Zeng_L6 <- read.delim(here("Data", "genelists", "Zeng.Layer.6.Markers.txt"), header=F) %>%
  as_tibble() %>%
  add_column(layer = 6)


Zeng_layer_markers <- rbind(Zeng_L1, Zeng_L2, Zeng_L3, Zeng_L4, Zeng_L5, Zeng_L6)
save(Zeng_layer_markers, file = here("Data", "genelists", "Zeng_layer_markers.Rdata"))
