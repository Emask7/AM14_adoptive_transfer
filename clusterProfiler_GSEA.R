# Make the gene list and sort the list in decreasing order ---------------------
  mk_list <- function(dat) {
    gene_list <- dat[, 2]
    names(gene_list) <- rownames(dat)
    gene_list <- sort(gene_list, decreasing = TRUE)
    return(gene_list)
  }
  
  mk_uniprot_list <- function(dat) {
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
  
  nrow(subset(res_PL23, padj <= 0.05))
  nrow(subset(res_PL23, padj <= 0.05 & abs(log2FoldChange) >= 0.6))
  
  
  geneList_PL23 <- mk_list(subset(res_PL23, padj <= 0.05))
  head(geneList_PL23)
  length(geneList_PL23)
  
  uni_geneList_PL23 <- mk_uniprot_list(subset(res_PL23, pvalue <= 0.05))
  head(uni_geneList_PL23)
  length(uni_geneList_PL23)
  

# Run GO gene set enrichment analysis ----------------------------------------- 
  go_PL23 <- gseGO(geneList = geneList_PL23, OrgDb = "org.Mm.eg.db", 
                   ont = "BP", keyType = "ENSEMBL", 
                   minGSSize = 10, maxGSSize = 500, eps = 1e-10, 
                   pvalueCutoff = 1, pAdjustMethod = "BH", verbose = TRUE)
  
  head(go_PL23)
  min(go_PL23$pvalue)
  min(go_PL23$p.adjust)
  max(go_PL23$p.adjust)
  
  # min(go_PL23$pvalue)
  # [1] 0.007846611
  # >   min(go_PL23$p.adjust)
  # [1] 0.5210378
  
  
  go_min3_PL23 <- gseGO(geneList = geneList_PL23, OrgDb = "org.Mm.eg.db", 
                   ont = "BP", keyType = "ENSEMBL", 
                   minGSSize = 3, maxGSSize = 500, eps = 1e-10, 
                   pvalueCutoff = 1, pAdjustMethod = "BH", verbose = TRUE)
  
  head(go_min3_PL23)
  min(go_min3_PL23$pvalue)
  min(go_min3_PL23$p.adjust)

  # >   min(go_min3_PL23$pvalue)
  # [1] 0.004006858
  # >   min(go_min3_PL23$p.adjust)
  # [1] 0.9210344
  
  
  go_max800_min3_PL23 <- gseGO(geneList = geneList_PL23, OrgDb = "org.Mm.eg.db", 
                        ont = "BP", keyType = "ENSEMBL", 
                        minGSSize = 3, maxGSSize = 800, eps = 1e-10, 
                        pvalueCutoff = 1, pAdjustMethod = "BH", verbose = TRUE)
  
  head(go_max800_min3_PL23)
  min(go_max800_min3_PL23$pvalue)
  min(go_max800_min3_PL23$p.adjust)
  
  
  
  
  go_KEGG_PL23 <- gseKEGG(geneList = uni_geneList_PL23, organism = "mmu", 
                          keyType = "uniprot", 
                          minGSSize = 10, maxGSSize = 500, eps = 1e-10, 
                          pvalueCutoff = 0.1, pAdjustMethod = "BH", verbose = TRUE)
  head(go_KEGG_PL23)
  min(go_KEGG_PL23$pvalue)
  min(go_KEGG_PL23$p.adjust)
  max(go_KEGG_PL23$p.adjust)
  
