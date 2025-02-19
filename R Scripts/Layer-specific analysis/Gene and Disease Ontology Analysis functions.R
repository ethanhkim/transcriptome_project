# Necessary libraries
library(cowplot)
library(readr)
library(ggplot2)
library(tidyr)
library(magrittr)
library(AnnotationDbi)
library(annotate)
library(GO.db)
library(tmod)
library(org.Hs.eg.db) 
library(here)
library(dplyr)
library(tibble)
library(metaseqR)
library(conflicted)
library(rvest)
library(httr)


# Set preferred packages
conflict_prefer("paste", "base")
conflict_prefer("desc", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("arrange", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("strsplit", "base")
conflict_prefer("content", "httr")


###############################################################################

# Functions to use in analysis ----

# Adapted from Leon French - code to create gene sets from tmod using DisGeNet Disease Ontology annotations
create_geneBackground_DO <- function(source_data) {
  geneBackground <- source_data %>%
    dplyr::select(gene_symbol) %>%
    pull ()
  names(geneBackground) <- geneBackground
  geneBackground <- sort(unique(names(as.list(geneBackground))))
}

# Create geneSetsTable for analysis using DisGeNet Disease Ontology annotations
geneSets_tmod_DO <- function(dataset_geneBackground_DO) {
  disgenet <- read_tsv(paste0(here("Data", "genelists", "DisGeNet", "curated_gene_disease_associations.tsv"))) 
  disgenet %<>% dplyr::select(symbol = geneSymbol, name = diseaseName, ID = diseaseId)
  disgenet %<>% filter(symbol %in% dataset_geneBackground_DO)
  geneLists <- group_by(disgenet, ID) %>% 
    dplyr::summarise(name = paste(unique(name), collapse = ","), genes = unique(list(symbol)), size = dplyr::n()) %>% 
    filter(size >= minGOgroupSize & size <= maxGOgroupSize) 
  namedLists <- geneLists$genes
  names(namedLists) <- geneLists$ID
  idToName <- data.frame(ID = geneLists$ID, Title = geneLists$name)
  geneSetsDO <- makeTmod(modules = idToName, modules2genes = namedLists)
  return(geneSetsDO)
}

# Create geneSetsTable for analysis for Gene Ontology annotations
geneSets_tmod_GO <- function(dataset_geneBackground) {
  if (exists("geneSetsGO") && length(geneSetsGO$MODULES2GENES) > 1000 ) { #assume it's already loaded - needs a fix to see if the variable is declared
  } else {
    go_object <- as.list(org.Hs.egGO2ALLEGS)
    symbolsInGO <- getSYMBOL(unique(unlist(go_object)), data=goSource)
    
    #build GO sets for tmod -slow
    tmodNames <- data.frame()
    modules2genes <- list()
    goGroupName <- names(go_object)[1]
    showMethods(Term)
    
    goCount <- length(go_object)
    count <- 1
    for(goGroupName in names(go_object)) {
      if (count %% 1000 == 0) print(paste(count, "of", goCount))
      count <- count + 1
      
      goGroup <- go_object[goGroupName]
      geneIDs <- unique(unlist(goGroup, use.names=F))  #discard evidence codes
      genesymbols <- unique(getSYMBOL(geneIDs, data=goSource))
      
      genesymbols <- intersect(genesymbols, dataset_geneBackground) #get size after intersecting with our full gene set
      if (!(length(genesymbols) >= minGOgroupSize & length(genesymbols) <= maxGOgroupSize)) next();
      
      modules2genes[goGroupName] <- list(genesymbols)
      tmodNames <- rbind(tmodNames, data.frame(ID=goGroupName, Title = Term(goGroupName)))
    }
    geneSetsGO <- makeTmod(modules = tmodNames, modules2genes = modules2genes)
  }
  return(geneSetsGO)
}

# Function to generate an ordered list of genes arranged by descending order in 
# terms of gene expression per layer
generate_AUC_genelist <- function(data, layer) {
  layer_label <- layer
  data %<>%
    select(gene_symbol:L6)
  sorted_data <- data[order(data[,layer_label], decreasing = TRUE),]
  genelist <- sorted_data %>% pull(gene_symbol)
  return(genelist)
}

#Function to create table of AUC and p-values associated with disease ontology terms
AUC_table <- function(AUC_genelist, geneSet) {
  #Creates a top and bottom 20 list of ontology definitions using AUC
  AUC_table <- tibble::as_tibble(tmodUtest(c(AUC_genelist), mset=geneSet, qval = 1.01, filter = T))
  AUC_table %<>% rowwise() %>% mutate(P.Value = P.Value * 2) %>% ungroup() %>% mutate(adj.P.Val=p.adjust(P.Value, method="fdr")) #tmod runs one-sided tests
  AUC_table %<>% rowwise() %>% mutate(aspect = Ontology(ID)) #add the source ontology (could be filtered for just biological process)
  
  #collapse genesets that have the exact same set of genes
  AUC_table %<>% rowwise() %>% mutate(genes = paste(sort(unlist(geneSet$MODULES2GENES[ID])), collapse = " "))
  AUC_table %<>% ungroup() %>% group_by(genes, N1) %>% arrange(Title) %>% 
    summarize(MainTitle = dplyr::first(Title),  ID=paste(ID, collapse=","), AUC = dplyr::first(AUC), 
              P.Value= dplyr::first(P.Value), aspect= dplyr::first(aspect), 
              otherNames = if_else(length(unique(Title)) > 1, paste(Title[2:length(Title)], collapse=", "), ""))
  AUC_table %<>% ungroup() %>% mutate(adj.P.Val=p.adjust(P.Value, method="fdr"))
  AUC_table %<>% dplyr::select(-genes)
  AUC_table %<>% arrange(P.Value)
  AUC_table$rank <- 1:nrow(AUC_table)
  AUC_table %<>% dplyr::select(MainTitle, N1, AUC, P.Value, adj.P.Val, everything()) %>% arrange(adj.P.Val)
}

# Widen scRNA-seq data by cell type
widen_scRNA_seq <- function(dataset) {
  wide_data <- dataset %>%
    group_by(class_label) %>%
    nest()
  return(wide_data)
}

# Function to generate a dataframe of sorted gene lists by layer for AUC
AUC_genelist_fn <- function(dataset) {
  # List of layers to create AUC genelists and tables for
  layer_list <- c("L1", "L2", "L3", "L4", "L5", "L6")
  
  AUC_genelist <- list()
  for (i in layer_list) {
    AUC_genelist[[i]] <- generate_AUC_genelist(dataset, i)
  }
  AUC_genelist_df <- as.data.frame(do.call(rbind, AUC_genelist)) %>%
    t() %>% as.data.frame()
}

# Function to generate list of tables of AUC values using sorted gene lists
AUC_table_list_fn <- function(AUC_genelist_df, geneSet) {
  AUC_table_list <- list()
  for (i in 1:6) {
    AUC_table_list[[i]] <- AUC_table(AUC_genelist_df[,i], geneSet)
  }
  return(AUC_table_list)
}

#Function to combine p-values across datasets joined by the same disease ontology term
combine_pval_AUC_table <- function(Maynard_table, He_table, Allen_table = NULL) {
  
  if (is.null(Allen_table)) {
    combined_table <- inner_join(Maynard_table, He_table, by = "MainTitle",
                                 suffix = c("_Maynard", "_He")) %>%
      column_to_rownames(var = "MainTitle") %>%
      dplyr::select(contains("P.value"), contains("AUC")) %>%
      filter(P.Value_He <= 1) %>%
      filter(P.Value_Maynard <= 1)
  } else {
    combined_table <- inner_join(Maynard_table, He_table, by = "MainTitle",
                                 suffix = c("_Maynard", "_He"))
    combined_table <- inner_join(combined_table, Allen_table, by = "MainTitle") %>%
      column_to_rownames(var = "MainTitle") %>%
      rename(P.Value_Allen = P.Value, AUC_Allen = AUC, N1_Allen = N1, adj.P.Val_Allen = adj.P.Val, ID_Allen = ID, aspect_Allen = aspect,
             otherNames_Allen = otherNames, rank_Allen = rank) %>%
      mutate(mean_AUC = (AUC_Maynard+AUC_He+AUC_Allen)/3) %>%
      dplyr::select(contains("P.value"), contains("AUC"), contains("ID")) %>%
      filter(P.Value_He <= 1) %>%
      filter(P.Value_Maynard <= 1) %>%
      filter(P.Value_Allen <= 1)
  }
  
  p_value_fisher <- fisher.method(combined_table %>% select(P.Value_Maynard:P.Value_Allen), method = c("fisher"))
  #adjust the p-value using fdr
  p_value_fisher$adj_p_fdr <- p.adjust(p_value_fisher$p.value, method = "fdr")
  
  combined_table$combined_p_value <- p_value_fisher$p.value
  combined_table$adj_combined_p_value <- p_value_fisher$adj_p_fdr
  
  combined_table <- combined_table %>%
    rownames_to_column(var = "Ontology_Term") %>%
    arrange(adj_combined_p_value)
  combined_table$rank <- 1:nrow(combined_table)
  return(combined_table)
}

# Combine the AUC tables using combine_pval_AUC_table function
combined_AUC_table_fn <- function(Maynard_AUC_table, He_AUC_table, Allen_AUC_table, layer_range = c(1:6)) {
  combined_AUC_table_list <- list()
  for (i in layer_range) {
    combined_AUC_table_list[[i]] <- combine_pval_AUC_table(Maynard_AUC_table[[i]], He_AUC_table[[i]], Allen_AUC_table[[i]])
  }
  return(combined_AUC_table_list)
}

# Function to take geneBackground, geneSets, genelist dataframe and create AUC value table
dataset_AUC_summary <- function(source_dataset, ontology, geneBackground) {
  geneBackground <- create_geneBackground_DO(source_dataset)
  if (ontology == "disease") {
    geneSets <- geneSets_tmod_DO(geneBackground)
  } else {
    geneSets <- geneSets_tmod_GO(geneBackground)
  }
  AUC_genelist_df <- AUC_genelist_fn(source_dataset)
  AUC_table_list <- AUC_table_list_fn(AUC_genelist_df, geneSets)
  return(AUC_table_list)
}

# Create DO table for top 20 and bottom 20 DO terms
extract_term_and_pvalue <- function(table, number_of_terms, layer_range = c(1:6)) {
  table_list_up <- list()
  table_list_down <- list()
  for (i in layer_range) {
    table_list_up[[i]] <- table[[i]] %>%
      filter(((AUC_Maynard > 0.5) & (AUC_He > 0.5) & (AUC_Allen > 0.5)) & adj_combined_p_value < 0.05 & P.Value_Allen < 0.05) %>%
      head(n = number_of_terms) %>%
      # Filter for adj but send raw 
      select(rank, Ontology_Term, mean_AUC, AUC_Maynard, AUC_He, AUC_Allen, ID_Maynard, combined_p_value) %>%
      rename(ID = ID_Maynard) %>%
      mutate(ID = gsub(",.*", "", ID)) %>%
      add_column(direction = "positive", layer = i) 
  }
  for (i in layer_range) {
    table_list_down[[i]] <- table[[i]] %>%
      filter(((AUC_Maynard < 0.5) & (AUC_He < 0.5) & (AUC_Allen < 0.5)) & adj_combined_p_value < 0.05 & P.Value_Allen < 0.05) %>%
      head(n = number_of_terms) %>%
      # Filter for adj but send raw 
      select(rank, Ontology_Term, mean_AUC, AUC_Maynard, AUC_He, AUC_Allen, ID_Maynard, combined_p_value) %>%
      rename(ID = ID_Maynard) %>%
      mutate(ID = gsub(",.*", "", ID)) %>%
      add_column(direction = "negative", layer = i) 
  }
  table_list_up <- do.call(rbind, table_list_up)
  table_list_down <- do.call(rbind, table_list_down)
  
  table_total_list <- rbind(table_list_up, table_list_down)
  return(table_total_list) 
}

# Create ReviGO table
extract_GO_and_pvalue <- function(table, number_of_terms, layer_range = c(1:6)) {
  
  reviGO_table_list_up <- list()
  reviGO_table_list_down <- list()
  for (i in layer_range) {
    reviGO_table_list_up[[i]] <- table[[i]] %>%
      filter((AUC_Maynard > 0.5) & (AUC_He > 0.5) & (AUC_Allen > 0.5) & adj_combined_p_value < 0.05 & P.Value_Allen < 0.05) %>%
      head(n = number_of_terms) %>%
      # Filter for adj but send raw 
      select(ID_Maynard, combined_p_value) %>%
      rename(ID = ID_Maynard) %>%
      mutate(ID = gsub(",.*", "", ID)) %>%
      add_column(direction = "positive", layer = i) 
  }
  for (i in layer_range) {
    reviGO_table_list_down[[i]] <- table[[i]] %>%
      filter((AUC_Maynard < 0.5) & (AUC_He < 0.5) & (AUC_Allen < 0.5) & adj_combined_p_value < 0.05 & P.Value_Allen < 0.05) %>%
      head(n = number_of_terms) %>%
      select(ID_Maynard, combined_p_value) %>%
      rename(ID = ID_Maynard) %>%
      mutate(ID = gsub(",.*", "", ID)) %>%
      add_column(direction = "negative", layer = i)
  }
  reviGO_table_list_up <- do.call(rbind, reviGO_table_list_up)
  reviGO_table_list_down <- do.call(rbind, reviGO_table_list_down)
  
  reviGO_table_total_list <- rbind(reviGO_table_list_up, reviGO_table_list_down)
  return(reviGO_table_total_list) 
}


run_revigo <- function(goID_and_pvalue_table) {
  Sys.sleep(.5) #delay half a second to prevent overloading revigo
  if(ncol(goID_and_pvalue_table) != 2 ) stop('Input data frame is not two columns wide, should be go ID and p-value only (eg. GO:0009268 1e-14)')
  url <- "http://revigo.irb.hr/revigo.jsp"
  fd <- list(
    submit ="submitToRevigo",
    goList  = format_tsv(goID_and_pvalue_table, col_names = F),
    #goList  = "GO:0006614	2.139547041149198e-33
    #GO:0006613	6.532814810756933e-33",
    #goList  = "GO:0006613 1e-14
    #GO:0006614 1e-14",
    cutoff = "0.70",
    isPValue = "yes",
    whatIsBetter = "higher",
    goSizes = "0", #0 representing "whole UniProt (default)", alternatives are "Homo sapiens"= 9606 or another NCBI species code
    measure = "SIMREL" 
  )
  h1 <- handle('') #reset the handle
  resp <- POST(url, body=fd, encode="form",  handle=h1) #use verbose() for checking
  
  #for debugging
  #html(content(resp, as = "text"))
  #write(content(resp, as = "text"), "/Users/lfrench/Downloads/z.txt")
  #write(content(resp, as = "text"), "/Users/lfrench/Downloads/z.html")
  
  #clean up of the table
  extracted_data_frame <- html_table(content(resp))[[1]]
  extracted_data_frame <- extracted_data_frame[-1,]
  colnames(extracted_data_frame) <- extracted_data_frame[1,]
  extracted_data_frame <- extracted_data_frame[-1,]
  extracted_data_frame %<>% as_tibble()
  extracted_data_frame %<>% mutate(frequency = gsub(" %", "", frequency))
  extracted_data_frame %<>% type.convert()
  extracted_data_frame
}

# From Leon French - code to run ReviGO on localized table
run_revigo <- function(goID_and_pvalue_table, layer_number, gene_direction) {
  
  goID_and_pvalue_table %<>%
    filter(layer == layer_number) %>%
    filter(direction == gene_direction) %>%
    select(ID, combined_p_value)
  
  Sys.sleep(.5) #delay half a second to prevent overloading revigo
  if(ncol(goID_and_pvalue_table) != 2 ) stop('Input data frame is not two columns wide, should be go ID and p-value only (eg. GO:0009268 1e-14)')
  url <- "http://revigo.irb.hr/revigo.jsp"
  fd <- list(
    submit ="submitToRevigo",
    goList  = format_tsv(goID_and_pvalue_table, col_names = F),
    #goList  = "GO:0006614	2.139547041149198e-33
    #GO:0006613	6.532814810756933e-33",
    #goList  = "GO:0006613 1e-14
    #GO:0006614 1e-14",
    cutoff = "0.70",
    isPValue = "yes",
    whatIsBetter = "higher",
    goSizes = "0", #0 representing "whole UniProt (default)", alternatives are "Homo sapiens"= 9606 or another NCBI species code
    measure = "SIMREL" 
  )
  h1 <- handle('') #reset the handle
  resp <- POST(url, body=fd, encode="form",  handle=h1) #use verbose() for checking
  
  #for debugging
  #html(content(resp, as = "text"))
  #write(content(resp, as = "text"), "/Users/lfrench/Downloads/z.txt")
  #write(content(resp, as = "text"), "/Users/lfrench/Downloads/z.html")
  
  #clean up of the table
  extracted_data_frame <- html_table(content(resp))[[1]]
  extracted_data_frame <- extracted_data_frame[-1,]
  colnames(extracted_data_frame) <- extracted_data_frame[1,]
  extracted_data_frame <- extracted_data_frame[-1,]
  extracted_data_frame %<>% as_tibble()
  extracted_data_frame %<>% mutate(frequency = gsub(" %", "", frequency))
  extracted_data_frame %<>% type.convert()
  extracted_data_frame
}

# Print common ontology terms ---
print_ontology_terms <- function(cell_specific_AUC_table, layer_number, direction) {
  if (direction == "positive") {
    print(head(filter(cell_specific_AUC_table[[layer_number]], 
                      (( (AUC_Maynard > 0.5) & (AUC_He > 0.5) & (AUC_Allen > 0.5)) & 
                         adj_combined_p_value < 0.05 & P.Value_Allen < 0.05)), n=20) %>%
            select(Ontology_Term:AUC_Allen, mean_AUC, combined_p_value, adj_combined_p_value,
                   rank))
  } else if (direction == "negative") {
    print(head(filter(cell_specific_AUC_table[[layer_number]], 
                      (( (AUC_Maynard < 0.5) & (AUC_He < 0.5) & (AUC_Allen < 0.5)) & 
                         adj_combined_p_value < 0.05 & P.Value_Allen < 0.05)), n=20),
          select(Ontology_Term:AUC_Allen, mean_AUC, combined_p_value, adj_combined_p_value,
                 rank))
  }
}


