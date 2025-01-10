# # Set the desired organism name here -------------------------------------------
#   organism = "org.Mm.eg.db"
#   BiocManager::install("org.Mm.eg.db")
#   library(organism, character.only = TRUE)
#   
# # Prepare Input
#   # reading in data from deseq2
#   df = read.csv("drosphila_example_de.csv", header=TRUE)
#   
#   
#   # PL23_filt_res <- data.frame(subset(res_PL23, !is.na(padj)))
#   # R848_filt_res <- data.frame(subset(res_R848, !is.na(padj)))
#   # PL23vR848_filt_res <- data.frame(subset(res_PL23vR848, !is.na(padj)))
#   
#   
#   # we want the log2 fold change 
#   original_gene_list <- df$log2FoldChange
#   
#   # name the vector
#   names(original_gene_list) <- df$X
# 
# # omit any NA values 
# gene_list<-na.omit(original_gene_list)
  
# Make the gene list and sort the list in decreasing order --------------------
  mk_list <- function(dat) {
    # df <- subset(dat, padj <= 0.05)
    gene_list <- dat[, 2]
    names(gene_list) <- rownames(dat)
    gene_list <- sort(gene_list, decreasing = TRUE)
    return(gene_list)
  }
  
  mk_uniprot_list <- function(dat) {
    # df <- subset(dat, padj <= 0.05)
    df <- data.frame(rownames(dat), dat[, 2])
    colnames(df) <- c("ENSEMBL", "log2FoldChange")
    uniprot_names <- bitr(df[, 1], fromType = "ENSEMBL", toType = "UNIPROT", 
                          OrgDb = org.Mm.eg.db, drop = FALSE)
    uni_df <- full_join(uniprot_names, df, by = "ENSEMBL")
    uni_df <- subset(uni_df, !is.na(UNIPROT))
    gene_list <- uni_df$log2FoldChange
    names(gene_list) <- uni_df$UNIPROT
    gene_list <- sort(gene_list, decreasing = TRUE)
    return(gene_list)
  }
  
  

  uni_df <- full_join(uniprot_names, df, by = "SYMBOL")
  head(uni_df)
  nrow(uni_df)
  
  geneList_PL23 <- mk_list(subset(res_PL23, padj <= 0.05))
  head(geneList_PL23)
  length(geneList_PL23)
  
  
  uni_geneList_PL23 <- mk_uniprot_list(subset(res_PL23, pvalue <= 0.01))
  head(uni_geneList_PL23)
  length(uni_geneList_PL23)
  
  nrow(subset(res_PL23, pvalue <= 0.01))
  
# Run GO gene set enrichment analysis ----------------------------------------- 
  go_PL23 <- gseGO(geneList = geneList_PL23, OrgDb = "org.Mm.eg.db", 
                   ont = "BP", keyType = "ENSEMBL", 
                   minGSSize = 10, maxGSSize = 500, eps = 1e-10, 
                   pvalueCutoff = 1, pAdjustMethod = "BH", verbose = TRUE)
  
  head(go_PL23)
  min(go_PL23$pvalue)
  min(go_PL23$p.adjust)
  max(go_PL23$p.adjust)
  
  
  
  go_KEGG_PL23 <- gseKEGG(geneList = uni_geneList_PL23, organism = "mmu", 
                          keyType = "uniprot", 
                          minGSSize = 10, maxGSSize = 500, eps = 1e-10, 
                          pvalueCutoff = 0.1, pAdjustMethod = "BH", verbose = TRUE)
  head(go_KEGG_PL23)
  min(go_KEGG_PL23$pvalue)
  min(go_KEGG_PL23$p.adjust)
  max(go_KEGG_PL23$p.adjust)
  
