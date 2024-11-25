# Set the desired organism name here -------------------------------------------
  organism = "org.Mm.eg.db"
  BiocManager::install(organism, character.only = TRUE)
  library(organism, character.only = TRUE)
  
# Prepare Input
  # reading in data from deseq2
  df = read.csv("drosphila_example_de.csv", header=TRUE)
  
  
  # PL23_filt_res <- data.frame(subset(res_PL23, !is.na(padj)))
  # R848_filt_res <- data.frame(subset(res_R848, !is.na(padj)))
  # PL23vR848_filt_res <- data.frame(subset(res_PL23vR848, !is.na(padj)))
  
  
  # we want the log2 fold change 
  original_gene_list <- df$log2FoldChange
  
  # name the vector
  names(original_gene_list) <- df$X
  
  # omit any NA values 
  gene_list<-na.omit(original_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  