#Maynard et al comparisons --
library(vegan)
library(dplyr)
library(tidyverse)

Maynard_dataset_unique <- Maynard_dataset %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  as_tibble()

## Separate by sample ID
Maynard_151507 <- Maynard_dataset_unique %>%
  dplyr::select("gene_symbol", contains("151507")) %>%
  column_to_rownames(var = "gene_symbol")
Maynard_151508 <- Maynard_dataset_unique %>%
  dplyr::select("gene_symbol", contains("151508")) %>%
  column_to_rownames(var = "gene_symbol")
Maynard_151509 <- Maynard_dataset_unique %>%
  dplyr::select("gene_symbol", contains("151509")) %>%
  column_to_rownames(var = "gene_symbol")
Maynard_151510 <- Maynard_dataset_unique %>%
  dplyr::select("gene_symbol", contains("151510")) %>%
  column_to_rownames(var = "gene_symbol")

Maynard_151673 <- Maynard_dataset_unique %>%
  dplyr::select("gene_symbol", contains("151673")) %>%
  column_to_rownames(var = "gene_symbol")
Maynard_151674 <- Maynard_dataset_unique %>%
  dplyr::select("gene_symbol", contains("151674")) %>%
  column_to_rownames(var = "gene_symbol")
Maynard_151675 <- Maynard_dataset_unique %>%
  dplyr::select("gene_symbol", contains("151675")) %>%
  column_to_rownames(var = "gene_symbol")
Maynard_151676 <- Maynard_dataset_unique %>%
  dplyr::select("gene_symbol", contains("151676")) %>%
  column_to_rownames(var = "gene_symbol")


Maynard_mantel <- function(x, y) {
  temp_dataset_1 <- x %>%
    t()
  temp_matrix_1 <- cor(temp_dataset_1, temp_dataset_1, method = "pearson")
  temp_matrix_1[is.na(temp_matrix_1)] = 0
  
  temp_dataset_2 <- y %>%
    t()
  temp_matrix_2 <- cor(temp_dataset_2, temp_dataset_2, method = "pearson")
  temp_matrix_2[is.na(temp_matrix_2)] = 0
  
  mantel(temp_dataset_1, temp_dataset_2, permutations = 999, method = "spearman") %>%
    print()
}

Maynard_507_508 <- Maynard_mantel(Maynard_151507, Maynard_151508) 
Maynard_507_509 <- Maynard_mantel(Maynard_151507, Maynard_151509)
Maynard_507_510 <- Maynard_mantel(Maynard_151507, Maynard_151510)
Maynard_507_673 <- Maynard_mantel(Maynard_151507, Maynard_151673)
Maynard_507_674 <- Maynard_mantel(Maynard_151507, Maynard_151674)
Maynard_507_675 <- Maynard_mantel(Maynard_151507, Maynard_151675)
Maynard_507_676 <- Maynard_mantel(Maynard_151507, Maynard_151676)

Maynard_508_509 <- Maynard_mantel(Maynard_151508, Maynard_151509) 
Maynard_508_510 <- Maynard_mantel(Maynard_151508, Maynard_151510)
Maynard_508_673 <- Maynard_mantel(Maynard_151508, Maynard_151673)
Maynard_508_674 <- Maynard_mantel(Maynard_151508, Maynard_151674)
Maynard_508_675 <- Maynard_mantel(Maynard_151508, Maynard_151675)
Maynard_508_676 <- Maynard_mantel(Maynard_151508, Maynard_151676)

Maynard_509_510 <- Maynard_mantel(Maynard_151509, Maynard_151510) 
Maynard_509_673 <- Maynard_mantel(Maynard_151509, Maynard_151673)
Maynard_509_674 <- Maynard_mantel(Maynard_151509, Maynard_151674)
Maynard_509_675 <- Maynard_mantel(Maynard_151509, Maynard_151675)
Maynard_509_676 <- Maynard_mantel(Maynard_151509, Maynard_151676)

Maynard_510_673 <- Maynard_mantel(Maynard_151510, Maynard_151673) 
Maynard_510_674 <- Maynard_mantel(Maynard_151510, Maynard_151674)
Maynard_510_675 <- Maynard_mantel(Maynard_151510, Maynard_151675)
Maynard_510_676 <- Maynard_mantel(Maynard_151510, Maynard_151676)

Maynard_673_674 <- Maynard_mantel(Maynard_151673, Maynard_151674) 
Maynard_673_675 <- Maynard_mantel(Maynard_151673, Maynard_151675)
Maynard_673_676 <- Maynard_mantel(Maynard_151673, Maynard_151676)

Maynard_674_675 <- Maynard_mantel(Maynard_151674, Maynard_151675) 
Maynard_674_676 <- Maynard_mantel(Maynard_151674, Maynard_151676)

Maynard_675_676 <- Maynard_mantel(Maynard_151675, Maynard_151676) 



Maynard_mantel_results <- tibble(
  Mantel_statistic = c(Maynard_507_508$statistic, Maynard_507_509$statistic, Maynard_507_510$statistic, Maynard_507_673$statistic, Maynard_507_674$statistic,
                      Maynard_507_675$statistic, Maynard_507_676$statistic, Maynard_508_509$statistic, Maynard_508_510$statistic, Maynard_508_673$statistic,
                      Maynard_508_674$statistic, Maynard_508_675$statistic, Maynard_508_676$statistic, Maynard_509_510$statistic, Maynard_509_673$statistic,
                      Maynard_509_674$statistic, Maynard_509_675$statistic, Maynard_509_676$statistic, Maynard_510_673$statistic, Maynard_510_674$statistic,
                      Maynard_510_675$statistic, Maynard_510_676$statistic, Maynard_673_674$statistic, Maynard_673_675$statistic, Maynard_673_676$statistic,
                      Maynard_674_675$statistic, Maynard_674_676$statistic, Maynard_675_676$statistic),
  Compared_datasets = c("507_508", "507_509", "507_510", "507_673", "507_674", "507_675", "507_676", "508_509", 
                        "508_510", "508_673", "508_674", "508_675", "508_676", "509_510", "509_673", "509_674",
                        "509_675", "509_676", "510_673", "510_674", "510_675", "510_676", "673_674", "673_675",
                        "673_676", "674_675", "674_676", "675_676")
)

hist(Maynard_mantel_results$Mantel_statistic)

ggplot(data = Maynard_mantel_results, aes(x = Compared_datasets, y = Mantel_statistic)) +
  geom_col() +
  theme_bw() +
  theme(panel.background = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45,
                                   margin = margin(b = 20),
                                   vjust = 0.25))



