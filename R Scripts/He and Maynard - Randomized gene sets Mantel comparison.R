## Randomize datasets ##

#Generate a random list of 931 genes, setting a different seed each time
#Use that list to extract those genes from He and Maynard, along with its values
#run cor for those values
#run mantel against base correlation matrix
#store correlation value in vector

#Function to create 1 random Mantel test statistic

generate_mantel_statistic <- function(x) {
  
  #Generate random list of 931 genes from the Maynard and He datasets based on random seed
  He_Maynard_common_gene_list <- intersect(He_DS1_Human$gene_symbol, Maynard_dataset$gene_symbol)
  set.seed(x, kind = NULL) 
  random_geneList <- sample(He_Maynard_common_gene_list, size = 931)
  
  #Generate datasets of layers 1 to 7 (WM) based on random list
  temp_He_dataset <-  He_DS1_Human %>%
    filter(gene_symbol %in% random_geneList) %>%
    distinct(gene_symbol, .keep_all = TRUE) %>%
    slice(match(random_geneList, gene_symbol)) %>%
    column_to_rownames(var = "gene_symbol") %>%
    t()
  
  temp_Maynard_dataset <- Maynard_dataset %>%
    as.data.frame(x) %>%
    filter(gene_symbol %in% random_geneList) %>%
    distinct(gene_symbol, .keep_all = TRUE) %>%
    slice(match(random_geneList, gene_symbol)) %>%
    column_to_rownames(var = "gene_symbol") %>%
    t()
  
  #Generate matrix
  temp_He_cor_matrix <- cor(temp_He_dataset, temp_He_dataset, method = "pearson")
  temp_Maynard_cor_matrix <- cor(temp_Maynard_dataset, temp_Maynard_dataset, method = "pearson")
  
  #Perform mantel test
  mantel_results <- mantel(temp_He_cor_matrix, temp_Maynard_cor_matrix, method = "spearman", permutations = 1, na.rm = TRUE)
  cat(mantel_results$statistic, "\n")
}


#Generate 100 random mantel test statistics from seeds 1:100
random_mantel_statistic_list <- as_tibble(capture.output(for (i in 1:100) {
  generate_mantel_statistic(i)
}))

random_mantel_statistic_list <- as.numeric(unlist(random_mantel_statistic_list)) %>%
  as_tibble() %>%
  add_column(seed = c(1:100)) 

random_mantel_statistic_list <- random_mantel_statistic_list %>%
  slice(-101)

hist(random_mantel_statistic_list$value)

ggplot(data = random_mantel_statistic_list, aes(x = seed, y = value)) +
  geom_bar(stat = "identity")
