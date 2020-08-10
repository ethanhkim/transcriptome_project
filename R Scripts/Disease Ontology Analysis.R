library(cowplot)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(magrittr)
library(AnnotationDbi)
library(annotate)
library(GO.db)
library(tmod)
library(org.Hs.eg.db) 
library(here)

#Input parameters:

species <- "human"
maxGOgroupSize <- 200
minGOgroupSize <- 10
methodToUse <- "AUC" #either "AUC" or "hyper" depending on AUC or hypergeometric test
datasetToUse <- "Maynard" #either "Maynard" or "He" depending on if analyzing Maynard et al, or He et al data, respectively
databaseToUse <- "DisGeNet" 
goSource <- 'org.Hs.eg'


# Adapted from Leon French - code to create gene sets from tmod using DisGeNet Disease Ontology annotations
disgenet <- read_tsv(paste0(here("Data", "genelists", "DisGeNet", "curated_gene_disease_associations.tsv"))) 
disgenet %<>% dplyr::select(symbol = geneSymbol, name = diseaseName, ID = diseaseId)
disgenet %<>% filter(symbol %in% geneBackground)
  
geneLists <- group_by(disgenet, ID) %>% 
  dplyr::summarise(name = paste(unique(name), collapse = ","), genes = unique(list(symbol)), size = dplyr::n()) %>% 
  filter(size >= minGOgroupSize & size <= maxGOgroupSize) 
namedLists <- geneLists$genes
names(namedLists) <- geneLists$ID
idToName <- data.frame(ID = geneLists$ID, Title = geneLists$name)
geneSets <- makeTmod(modules = idToName, modules2genes = namedLists)
geneSetsTable <- tmod2DataFrame(geneSets)
  

#Specifying layers is done within the script - find "Change per layer to examine:"
if (methodToUse == "AUC") {
  if (datasetToUse == "Maynard") {
    Maynard_genelist <- Maynard_dataset_average %>%
      arrange(desc(Layer_2)) %>% # Change per layer to examine: Layer_1, Layer_2, Layer_3, Layer_4, Layer_5, Layer_6
      dplyr::select(gene_symbol) %>%
      pull()
    names(Maynard_genelist) <- Maynard_genelist 
    sortedGenes <- unique(names(as.list(Maynard_genelist)))
  } else { #if datasetToUse == "He"
    He_genelist <- He_DS1_Human_averaged %>%
      arrange(desc(Layer_4)) %>% # Change per layer to examine: Layer_1, Layer_2, Layer_3, Layer_4, Layer_5, Layer_6
      dplyr::select(gene_symbol) %>%
      pull()
    names(He_genelist) <- He_genelist 
    sortedGenes <- unique(names(as.list(He_genelist)))
  }
  
  #Creates a top and bottom 20 list of ontology definitions using AUC
  result <- tbl_df(tmodUtest(c(sortedGenes), mset=geneSets, qval = 1.01, filter = T))
  result %<>% rowwise() %>% mutate(P.Value = P.Value * 2) %>% ungroup() %>% mutate(adj.P.Val=p.adjust(P.Value, method="fdr")) #tmod runs one-sided tests
  result %<>% rowwise() %>% mutate(aspect = Ontology(ID)) #add the source ontology (could be filtered for just biological process)
  
  #collapse genesets that have the exact same set of genes
  result %<>% rowwise() %>% mutate(genes = paste(sort(unlist(geneSets$MODULES2GENES[ID])), collapse = " "))
  result %<>% ungroup() %>% group_by(genes, N1) %>% arrange(Title) %>% 
    summarize(MainTitle = dplyr::first(Title),  ID=paste(ID, collapse=","), AUC = dplyr::first(AUC), 
              P.Value= dplyr::first(P.Value), aspect= dplyr::first(aspect), 
              otherNames = if_else(length(unique(Title)) > 1, paste(Title[2:length(Title)], collapse=", "), ""))
  result %<>% ungroup() %>% mutate(adj.P.Val=p.adjust(P.Value, method="fdr"))
  result %<>% dplyr::select(-genes)
  
  result %<>% arrange(P.Value)
  result$rank <- 1:nrow(result)
  result %<>% dplyr::select(MainTitle, N1, AUC, P.Value, adj.P.Val, everything()) %>% arrange(adj.P.Val)
  
  #print top 20 at top and bottom
  print(head(filter(result, AUC > 0.5) %>% dplyr::select(-ID), n=20))
  print(head(filter(result, AUC < 0.5) %>% dplyr::select(-ID), n=20))
  
  source("./ROCPlots.R")
  
  top_go_term <- result %>%
    filter(AUC > 0.5) %>% dplyr::select(ID) %>% slice_head() %>%
    mutate(ID = gsub(",.*","",ID)) %>% as.character() 
  
  bottom_go_term <- result %>%
    filter(AUC < 0.5) %>% dplyr::select(ID) %>% slice_head() %>%
    mutate(ID = gsub(",.*","",ID)) %>% as.character() 
  
  #filter sorted genes for those in GO
  sortedGenes <- intersect(sortedGenes, geneSets$GENES$ID)
  
  #inspect GO term and list what genes are included:
  geneSets$MODULES2GENES[bottom_go_term] %>% unlist %>% unname
  
  #plot the top and bottom GO Groups
  plots <- createPlots(sortedGenes, c(bottom_go_term, top_go_term), geneSets)
  plots$AUCPlot
  plots$rasterPlot
  (bothPlots <- plot_grid(plots$AUCPlot, plots$rasterPlot,  nrow = 2, align = "v", rel_heights=c(1,0.2),scale = 0.95))
  
} else { #if method == "hyper"
  if (datasetToUse == "Maynard") {
    Maynard_genelist <- Maynard_dataset_average %>%
      dplyr::select(gene_symbol, Layer_1_marker) %>%
      filter(!is.na(Layer_1_marker)) %>%
      dplyr::select(gene_symbol) %>%
      pull()
    names(Maynard_genelist) <- Maynard_genelist 
    
    hitListGenes <- sort(unique(names(as.list(Maynard_genelist))))
    
    sortedGenes <- Maynard_dataset_average %>%
      dplyr::select(gene_symbol) %>%
      pull()
    names(sortedGenes) <- sortedGenes
    sortedGenes <- sort(unique(names(as.list(sortedGenes))))
  } else { #if datasetToUse == "He"
    He_genelist <- He_DS1_Human_averaged %>%
      dplyr::select(gene_symbol, Layer_1_marker) %>%
      filter(!is.na(Layer_1_marker)) %>%
      dplyr::select(gene_symbol) %>%
      pull()
    names(He_genelist) <- He_genelist 
    
    hitListGenes <- sort(unique(names(as.list(He_genelist))))
    
    sortedGenes <- He_DS1_Human_averaged %>%
      dplyr::select(gene_symbol) %>%
      pull()
    names(sortedGenes) <- sortedGenes
    sortedGenes <- sort(unique(names(as.list(sortedGenes))))
  }
  
  
  #Create different geneBackground depending on the dataset specified
  if (datasetToUse == "Maynard") {
    geneBackground <- Maynard_dataset_average %>%
      dplyr::select(gene_symbol) %>%
      pull ()
    names(geneBackground) <- geneBackground
    geneBackground <- sort(unique(names(as.list(geneBackground))))
  } else {
    geneBackground <- He_DS1_Human_averaged %>%
      dplyr::select(gene_symbol) %>%
      pull ()
    names(geneBackground) <- geneBackground
    geneBackground <- sort(unique(names(as.list(geneBackground))))
  }
  
  
  #this just compares two lists and creates a top 20 list using a hypergeometric test
  result <- tbl_df(tmodHGtest(fg = hitListGenes, bg = sortedGenes, mset=geneSets, qval = 1.01, filter = T))
  
  result %<>% rowwise() %>% mutate(aspect = Ontology(ID)) #add the source ontology (could be filterd for just biological process)
  result$rank <- 1:nrow(result)
  result %<>% dplyr::select(Title, overlap = b, setSize=B, hitListSize = n, genesWithGOAnnotation = N, P.Value, adj.P.Val, everything()) %>% arrange(adj.P.Val)
  
  #print top 20 
  print(head(result %>% dplyr::select(-ID, -E, -rank), n=20))
}
