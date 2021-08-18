## Compare transposed Maynard datasets ##

library(dplyr)
library(tibble)
library(magrittr)
library(data.table)
library(conflicted) # Easily manage conflicting libraries
library(vegan) # Package for mantel test
library(here)
library(parallel) # Parallelize Mantel and WGCNA corr()
library(WGCNA) # Faster corr() than base
library(edgeR)

# Set conflicts
conflict_prefer('intersect', 'dplyr')
conflict_prefer("filter", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("cor", "WGCNA")

# Source mantel test function
source(here("R Scripts", "Mantel testing", "mantel_test.R"))

# Normalize per sample
normalize_sample <- function(sample_df, CPM_filter = FALSE) {
  
  df <- sample_df
  genes <- df$gene_symbol
  
  if (CPM_filter == TRUE) {
    logCPM_df <- df %>%
      select(-gene_symbol) %>%
      mutate_all(~. +1) %>%
      cpm() %>%
      as.data.frame() %>%
      add_column(gene_symbol = genes) %>%
      # Filter out samples of CPM < 0.1
      filter_at(vars(-gene_symbol), all_vars(. > .1)) %>%
      column_to_rownames(var = "gene_symbol")
    names(logCPM_df) <- gsub(x = names(logCPM_df), pattern = ".*_", replacement = "") 
    names(logCPM_df) <- gsub(x = names(logCPM_df), pattern = "Layer", replacement = "L")
    filtered_genes <- rownames(logCPM_df)
    logCPM_df %<>%
      map_df(log2) %>%
      add_column(gene_symbol = filtered_genes) %>%
      select(gene_symbol, L1, L2, L3, L4, L5, L6, WM) %>%
      select(gene_symbol, everything())
  } else {
    # Normalize data
    logCPM_df <- df %>%
      # Remove gene_symbol column for cpm()
      select(-gene_symbol) %>% 
      # Add +1 to remove 0's for log transformation
      mutate_all(~. +1) %>%
      # CPM normalize with log = T
      cpm(log = T) %>% as.data.frame()
    names(logCPM_df) <- gsub(x = names(logCPM_df), pattern = ".*_", replacement = "") 
    names(logCPM_df) <- gsub(x = names(logCPM_df), pattern = "Layer", replacement = "L")
    genes <- rownames(logCPM_df)
    logCPM_df %<>%
      # Add back in gene symbols
      add_column(gene_symbol = genes) %>%
      select(gene_symbol, everything())
  }
  
  return(logCPM_df)
}


## Load Maynard data ##
Maynard_layer_data <- fread(here("Data", "raw_data", "Maynard et al",
                                 "layer_level_data.csv")) %>%
  select(-V1)

# Zeng markers
load(here("Data", "processed_data", "Zeng_dataset_long.Rdata"))
Zeng_marker_genes <- unique(Zeng_dataset_long$gene_symbol)

# Filter Maynard by list of genes of CPM > 0.1 in bulked data

# Create list of tibbles for donor 1, separated by sample
Maynard_donor_1 <- list()
for (sample in c("151507", "151508", "151509", "151510")) {
  Maynard_donor_1[[sample]] <- Maynard_layer_data %>%
    select("gene_symbol", contains(sample)) %>%
    normalize_sample()
}

Maynard_donor_1_CPM <- list()
for (sample in c("151507", "151508", "151509", "151510")) {
  Maynard_donor_1_CPM[[sample]] <- Maynard_layer_data %>%
    select("gene_symbol", contains(sample)) %>%
    normalize_sample(CPM_filter = TRUE)
}



# Create list of tibbles for donor 2, separated by sample
Maynard_donor_2 <- list()
for (sample in c("151673", "151674", "151675", "151676")) {
  Maynard_donor_2[[sample]] <- Maynard_layer_data %>%
    select("gene_symbol", contains(sample)) %>%
    normalize_sample()
}


# Create list of tibbles for donor 2, separated by sample
Maynard_donor_2_CPM <- list()
for (sample in c("151673", "151674", "151675", "151676")) {
  Maynard_donor_2_CPM[[sample]] <- Maynard_layer_data %>%
    select("gene_symbol", contains(sample)) %>%
    normalize_sample(CPM_filter = TRUE)
}

CPM_filter <- TRUE

# Intra-donor 1 Mantel tests
if (CPM_filter == FALSE) {
  Maynard_507_508 <- mantel_test(Maynard_donor_1$'151507', 
                                 Maynard_donor_1$'151508') 
  Maynard_507_509 <- mantel_test(Maynard_donor_1$'151507', 
                                 Maynard_donor_1$'151509')
  Maynard_507_510 <- mantel_test(Maynard_donor_1$'151507', 
                                 Maynard_donor_1$'151510')
  Maynard_508_509 <- mantel_test(Maynard_donor_1$'151508', 
                                 Maynard_donor_1$'151509')
  Maynard_508_510 <- mantel_test(Maynard_donor_1$'151508', 
                                 Maynard_donor_1$'151510')
  Maynard_509_510 <- mantel_test(Maynard_donor_1$'151509', 
                                 Maynard_donor_1$'151510')
  
  donor_1_mantel <- tibble(
    comparison = c("507_508", "507_509", "507_510", "508_509",
                   "508_510", "509_510"),
    statistic = c(Maynard_507_508$statistic, Maynard_507_509$statistic,
                  Maynard_507_510$statistic, Maynard_508_509$statistic,
                  Maynard_508_510$statistic, Maynard_509_510$statistic),
    p_val = c(Maynard_507_508$signif, Maynard_507_509$signif,
              Maynard_507_510$signif, Maynard_508_509$signif,
              Maynard_508_510$signif, Maynard_509_510$signif)
  )
} else {
  Maynard_507_508 <- mantel_test(Maynard_donor_1_CPM$'151507', 
                                 Maynard_donor_1_CPM$'151508') 
  Maynard_507_509 <- mantel_test(Maynard_donor_1_CPM$'151507', 
                                 Maynard_donor_1_CPM$'151509')
  Maynard_507_510 <- mantel_test(Maynard_donor_1_CPM$'151507', 
                                 Maynard_donor_1_CPM$'151510')
  Maynard_508_509 <- mantel_test(Maynard_donor_1_CPM$'151508', 
                                 Maynard_donor_1_CPM$'151509')
  Maynard_508_510 <- mantel_test(Maynard_donor_1_CPM$'151508', 
                                 Maynard_donor_1_CPM$'151510')
  Maynard_509_510 <- mantel_test(Maynard_donor_1_CPM$'151509', 
                                 Maynard_donor_1_CPM$'151510')
  
  donor_1_mantel <- tibble(
    comparison = c("507_508", "507_509", "507_510", "508_509",
                   "508_510", "509_510"),
    statistic = c(Maynard_507_508$statistic, Maynard_507_509$statistic,
                  Maynard_507_510$statistic, Maynard_508_509$statistic,
                  Maynard_508_510$statistic, Maynard_509_510$statistic),
    p_val = c(Maynard_507_508$signif, Maynard_507_509$signif,
              Maynard_507_510$signif, Maynard_508_509$signif,
              Maynard_508_510$signif, Maynard_509_510$signif)
  )
}

if (CPM_filter == FALSE) {
  
  # Intra-donor 2 Mantel tests
  Maynard_673_674 <- mantel_test(Maynard_donor_2$'151673', 
                                 Maynard_donor_2$'151674')
  Maynard_673_675 <- mantel_test(Maynard_donor_2$'151673', 
                                 Maynard_donor_2$'151675')
  Maynard_673_676 <- mantel_test(Maynard_donor_2$'151673', 
                                 Maynard_donor_2$'151676')
  Maynard_674_675 <- mantel_test(Maynard_donor_2$'151674', 
                                 Maynard_donor_2$'151675')
  Maynard_674_676 <- mantel_test(Maynard_donor_2$'151674', 
                                 Maynard_donor_2$'151676')
  Maynard_675_676 <- mantel_test(Maynard_donor_2$'151675', 
                                 Maynard_donor_2$'151676')
  
  donor_2_mantel <- tibble(
    comparison = c("673_674", "673_675", "673_676", 
                   "674_675", "674_676", "675_676"),
    statistic = c(Maynard_673_674$statistic, Maynard_673_675$statistic,
                  Maynard_673_676$statistic, Maynard_674_675$statistic,
                  Maynard_674_676$statistic, Maynard_675_676$statistic)
  )
  
} else {
  
  # Intra-donor 2 Mantel tests
  Maynard_673_674 <- mantel_test(Maynard_donor_2_CPM$'151673', 
                                 Maynard_donor_2_CPM$'151674')
  Maynard_673_675 <- mantel_test(Maynard_donor_2_CPM$'151673', 
                                 Maynard_donor_2_CPM$'151675')
  Maynard_673_676 <- mantel_test(Maynard_donor_2_CPM$'151673', 
                                 Maynard_donor_2_CPM$'151676')
  Maynard_674_675 <- mantel_test(Maynard_donor_2_CPM$'151674', 
                                 Maynard_donor_2_CPM$'151675')
  Maynard_674_676 <- mantel_test(Maynard_donor_2_CPM$'151674', 
                                 Maynard_donor_2_CPM$'151676')
  Maynard_675_676 <- mantel_test(Maynard_donor_2_CPM$'151675', 
                                 Maynard_donor_2_CPM$'151676')
  
  donor_2_mantel <- tibble(
    comparison = c("673_674", "673_675", "673_676", 
                   "674_675", "674_676", "675_676"),
    statistic = c(Maynard_673_674$statistic, Maynard_673_675$statistic,
                  Maynard_673_676$statistic, Maynard_674_675$statistic,
                  Maynard_674_676$statistic, Maynard_675_676$statistic)
  )
  
}



## Inter-donor testing
# Sample 151507 vs. samples from donor 2
donor_1_507_comparison <- list()
for (sample in c("151673", "151674", "151675", "151676")) {
  donor_1_507_comparison[[sample]] <- mantel_test(Maynard_donor_1$'151507',
                                                  Maynard_donor_2$sample)
}
Maynard_507_673 <- mantel_test(Maynard_donor_1$'151507', Maynard_donor_2$'151673')
Maynard_507_674 <- mantel_test(Maynard_donor_1$'151507', Maynard_donor_2$'151674')
Maynard_507_675 <- mantel_test(Maynard_donor_1$'151507', Maynard_donor_2$'151675')
Maynard_507_676 <- mantel_test(Maynard_donor_1$'151507', Maynard_donor_2$'151676')

Maynard_508_673 <- mantel_test(Maynard_donor_1$'151508', Maynard_donor_2$'151673')
Maynard_508_674 <- mantel_test(Maynard_donor_1$'151508', Maynard_donor_2$'151674')
Maynard_508_675 <- mantel_test(Maynard_donor_1$'151508', Maynard_donor_2$'151675')
Maynard_508_676 <- mantel_test(Maynard_donor_1$'151508', Maynard_donor_2$'151676')

Maynard_509_673 <- mantel_test(Maynard_donor_1$'151509', Maynard_donor_2$'151673')
Maynard_509_674 <- mantel_test(Maynard_donor_1$'151509', Maynard_donor_2$'151674')
Maynard_509_675 <- mantel_test(Maynard_donor_1$'151509', Maynard_donor_2$'151675')
Maynard_509_676 <- mantel_test(Maynard_donor_1$'151509', Maynard_donor_2$'151676')

Maynard_510_673 <- mantel_test(Maynard_donor_1$'151510', Maynard_donor_2$'151673') 
Maynard_510_674 <- mantel_test(Maynard_donor_1$'151510', Maynard_donor_2$'151674')
Maynard_510_675 <- mantel_test(Maynard_donor_1$'151510', Maynard_donor_2$'151675')
Maynard_510_676 <- mantel_test(Maynard_donor_1$'151510', Maynard_donor_2$'151676')


intra_donor_results <- tibble(
  comparison = c("507_673", "507_674", "507_675", "507_676", 
                 "508_673", "508_674", "508_675", "508_676", 
                 "509_673", "509_674", "509_675", "509_676", 
                 "510_673", "510_674", "510_675", "510_676"),
  statistic = c(Maynard_507_673$statistic, Maynard_507_674$statistic, Maynard_507_675$statistic, Maynard_507_676$statistic, 
                Maynard_508_673$statistic, Maynard_508_674$statistic, Maynard_508_675$statistic, Maynard_508_676$statistic, 
                Maynard_509_673$statistic, Maynard_509_674$statistic, Maynard_509_675$statistic, Maynard_509_676$statistic, 
                Maynard_510_673$statistic, Maynard_510_674$statistic, Maynard_510_675$statistic, Maynard_510_676$statistic)
)


  

# Create list of tibbles for donor 1, separated by sample
Maynard_donor_1_Zeng <- list()
for (sample in c("151507", "151508", "151509", "151510")) {
  Maynard_donor_1_Zeng[[sample]] <- Maynard_layer_data %>%
    select("gene_symbol", contains(sample)) %>%
    normalize_sample(CPM_filter = TRUE)
}


# Create list of tibbles for donor 2, separated by sample
Maynard_donor_2_Zeng <- list()
for (sample in c("151673", "151674", "151675", "151676")) {
  Maynard_donor_2_Zeng[[sample]] <- Maynard_layer_data %>%
    select("gene_symbol", contains(sample)) %>%
    normalize_sample(CPM_filter = TRUE)
}

# Intra-donor 1 Mantel tests
Maynard_507_508 <- mantel_test(Maynard_donor_1_Zeng$'151507', 
                               Maynard_donor_1_Zeng$'151508') 
Maynard_507_509 <- mantel_test(Maynard_donor_1_Zeng$'151507', 
                               Maynard_donor_1_Zeng$'151509')
Maynard_507_510 <- mantel_test(Maynard_donor_1_Zeng$'151507', 
                               Maynard_donor_1_Zeng$'151510')
Maynard_508_509 <- mantel_test(Maynard_donor_1_Zeng$'151508', 
                               Maynard_donor_1_Zeng$'151509')
Maynard_508_510 <- mantel_test(Maynard_donor_1_Zeng$'151508', 
                               Maynard_donor_1_Zeng$'151510')
Maynard_509_510 <- mantel_test(Maynard_donor_1_Zeng$'151509', 
                               Maynard_donor_1_Zeng$'151510')

donor_1_mantel <- tibble(
  comparison = c("507_508", "507_509", "507_510", "508_509",
                 "508_510", "509_510"),
  statistic = c(Maynard_507_508$statistic, Maynard_507_509$statistic,
                Maynard_507_510$statistic, Maynard_508_509$statistic,
                Maynard_508_510$statistic, Maynard_509_510$statistic),
  p_val = c(Maynard_507_508$signif, Maynard_507_509$signif,
            Maynard_507_510$signif, Maynard_508_509$signif,
            Maynard_508_510$signif, Maynard_509_510$signif)
)

# Intra-donor 2 Mantel tests
Maynard_673_674 <- mantel_test(Maynard_donor_2_Zeng$'151673', 
                               Maynard_donor_2_Zeng$'151674')
Maynard_673_675 <- mantel_test(Maynard_donor_2_Zeng$'151673', 
                               Maynard_donor_2_Zeng$'151675')
Maynard_673_676 <- mantel_test(Maynard_donor_2_Zeng$'151673', 
                               Maynard_donor_2_Zeng$'151676')
Maynard_674_675 <- mantel_test(Maynard_donor_2_Zeng$'151674', 
                               Maynard_donor_2_Zeng$'151675')
Maynard_674_676 <- mantel_test(Maynard_donor_2_Zeng$'151674', 
                               Maynard_donor_2_Zeng$'151676')
Maynard_675_676 <- mantel_test(Maynard_donor_2_Zeng$'151675', 
                               Maynard_donor_2_Zeng$'151676')

donor_2_mantel <- tibble(
  comparison = c("673_674", "673_675", "673_676", 
                 "674_675", "674_676", "675_676"),
  statistic = c(Maynard_673_674$statistic, Maynard_673_675$statistic,
                Maynard_673_676$statistic, Maynard_674_675$statistic,
                Maynard_674_676$statistic, Maynard_675_676$statistic)
)


Maynard_507_673 <- mantel_test(Maynard_donor_1_Zeng$'151507', Maynard_donor_2_Zeng$'151673')
Maynard_507_674 <- mantel_test(Maynard_donor_1_Zeng$'151507', Maynard_donor_2_Zeng$'151674')
Maynard_507_675 <- mantel_test(Maynard_donor_1_Zeng$'151507', Maynard_donor_2_Zeng$'151675')
Maynard_507_676 <- mantel_test(Maynard_donor_1_Zeng$'151507', Maynard_donor_2_Zeng$'151676')

Maynard_508_673 <- mantel_test(Maynard_donor_1_Zeng$'151508', Maynard_donor_2_Zeng$'151673')
Maynard_508_674 <- mantel_test(Maynard_donor_1_Zeng$'151508', Maynard_donor_2_Zeng$'151674')
Maynard_508_675 <- mantel_test(Maynard_donor_1_Zeng$'151508', Maynard_donor_2_Zeng$'151675')
Maynard_508_676 <- mantel_test(Maynard_donor_1_Zeng$'151508', Maynard_donor_2_Zeng$'151676')

Maynard_509_673 <- mantel_test(Maynard_donor_1_Zeng$'151509', Maynard_donor_2_Zeng$'151673')
Maynard_509_674 <- mantel_test(Maynard_donor_1_Zeng$'151509', Maynard_donor_2_Zeng$'151674')
Maynard_509_675 <- mantel_test(Maynard_donor_1_Zeng$'151509', Maynard_donor_2_Zeng$'151675')
Maynard_509_676 <- mantel_test(Maynard_donor_1_Zeng$'151509', Maynard_donor_2_Zeng$'151676')

Maynard_510_673 <- mantel_test(Maynard_donor_1_Zeng$'151510', Maynard_donor_2_Zeng$'151673') 
Maynard_510_674 <- mantel_test(Maynard_donor_1_Zeng$'151510', Maynard_donor_2_Zeng$'151674')
Maynard_510_675 <- mantel_test(Maynard_donor_1_Zeng$'151510', Maynard_donor_2_Zeng$'151675')
Maynard_510_676 <- mantel_test(Maynard_donor_1_Zeng$'151510', Maynard_donor_2_Zeng$'151676')


intra_donor_results <- tibble(
  comparison = c("507_673", "507_674", "507_675", "507_676", 
                 "508_673", "508_674", "508_675", "508_676", 
                 "509_673", "509_674", "509_675", "509_676", 
                 "510_673", "510_674", "510_675", "510_676"),
  statistic = c(Maynard_507_673$statistic, Maynard_507_674$statistic, Maynard_507_675$statistic, Maynard_507_676$statistic, 
                Maynard_508_673$statistic, Maynard_508_674$statistic, Maynard_508_675$statistic, Maynard_508_676$statistic, 
                Maynard_509_673$statistic, Maynard_509_674$statistic, Maynard_509_675$statistic, Maynard_509_676$statistic, 
                Maynard_510_673$statistic, Maynard_510_674$statistic, Maynard_510_675$statistic, Maynard_510_676$statistic)
)