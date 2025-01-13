# Use for testing different ID types -------------------------------------------
  # PL23_symbols <- res_PL23
  # PL23_ensembl <- res_PL23
  PL23_entrez <- res_PL23
  
  summary_v2(PL23_symbols)
  summary_v2(PL23_ensembl)
  summary_v2(PL23_entrez)
  
  PL23vR848_entrez <- res_PL23vR848

# Make the gene list and sort the list in decreasing order ---------------------
  mk_list <- function(dat) {
    gene_list <- dat[, 2]
    names(gene_list) <- rownames(dat)
    return(sort(gene_list, decreasing = TRUE))
  }

  # nrow(subset(res_PL23, padj <= 0.05))
  # nrow(subset(res_PL23, padj <= 0.05 & abs(log2FoldChange) >= 0.6))
  # geneList_PL23 <- mk_list(subset(res_PL23, padj <= 0.05))
  # head(geneList_PL23)
  # length(geneList_PL23)

# Run GO gene set enrichment analysis ------------------------------------------ 
  symbolList_PL23 <- mk_list(subset(PL23_symbols, padj <= 0.05))
  go_symbols_PL23 <- gseGO(geneList = symbolList_PL23, keyType = "SYMBOL", 
                           OrgDb = "org.Mm.eg.db", ont = "BP", 
                           minGSSize = 10, maxGSSize = 500, 
                           pvalueCutoff = 1, pAdjustMethod = "BH")
  min(go_symbols_PL23$pvalue)
  min(go_symbols_PL23$p.adjust)
  nrow(subset(go_symbols_PL23, go_symbols_PL23$pvalue < 0.05))
  nrow(subset(go_symbols_PL23, go_symbols_PL23$p.adjust < 0.05))
  head(go_symbols_PL23)
  
  ensemblList_PL23 <- mk_list(subset(PL23_ensembl, padj <= 0.05))
  go_ensembl_PL23 <- gseGO(geneList = ensemblList_PL23, keyType = "ENSEMBL", 
                           OrgDb = "org.Mm.eg.db", ont = "BP", 
                           minGSSize = 10, maxGSSize = 500, 
                           pvalueCutoff = 1, pAdjustMethod = "BH")
  min(go_ensembl_PL23$pvalue)
  min(go_ensembl_PL23$p.adjust)
  nrow(subset(go_ensembl_PL23, go_ensembl_PL23$pvalue < 0.05))
  nrow(subset(go_ensembl_PL23, go_ensembl_PL23$p.adjust < 0.05))
  head(go_ensembl_PL23)
  
  entrezList_PL23 <- mk_list(subset(PL23_entrez, padj <= 0.05))
  go_entrez_PL23 <- gseGO(geneList = entrezList_PL23, keyType = "ENTREZID", 
                          OrgDb = "org.Mm.eg.db", ont = "BP", 
                          # minGSSize = 10, maxGSSize = 500,
                          minGSSize = 5, maxGSSize = 100,
                          pvalueCutoff = 1, pAdjustMethod = "BH")
  min(go_entrez_PL23$pvalue)
  min(go_entrez_PL23$p.adjust)
  nrow(go_entrez_PL23)
  nrow(subset(go_entrez_PL23, go_entrez_PL23$pvalue < 0.05))
  nrow(subset(go_entrez_PL23, go_entrez_PL23$p.adjust < 0.05))
  head(go_entrez_PL23)
  

# KEGG pathway over-representation analysis ----------------------------------
  enrichList_PL23 <- mk_list(subset(PL23_entrez, 
                                    padj <= 0.05 & abs(log2FoldChange) >= 0.6))
  enrichKEGG_PL23 <- enrichKEGG(gene = names(enrichList_PL23), 
                                organism = "mmu", keyType = "ncbi-geneid", 
                                # minGSSize = 10, maxGSSize = 500,
                                minGSSize = 5, maxGSSize = 100,
                                qvalueCutoff = 1,
                                pvalueCutoff = 1, pAdjustMethod = "BH")

  min(enrichKEGG_PL23$pvalue)
  min(enrichKEGG_PL23$p.adjust)
  nrow(subset(enrichKEGG_PL23, enrichKEGG_PL23$pvalue < 0.05))
  nrow(subset(enrichKEGG_PL23, enrichKEGG_PL23$p.adjust < 0.05))
  head(enrichKEGG_PL23)
  subset(enrichKEGG_PL23, enrichKEGG_PL23$p.adjust < 0.05)$Description
  subset(enrichKEGG_PL23, enrichKEGG_PL23$pvalue < 0.05)$Description
  
  subset(enrichKEGG_PL23, enrichKEGG_PL23$Description == "Systemic lupus erythematosus - Mus musculus (house mouse)")
  subset(enrichKEGG_PL23, enrichKEGG_PL23$Description == "NF-kappa B signaling pathway - Mus musculus (house mouse)")

  
  
  enrichList_PL23vR848 <- mk_list(subset(PL23vR848_entrez, 
                                         padj <= 0.05 & abs(log2FoldChange) >= 0.6))
  enrichKEGG_PL23vR848 <- enrichKEGG(gene = names(enrichList_PL23vR848), 
                                     organism = "mmu", keyType = "ncbi-geneid", 
                                     minGSSize = 5, maxGSSize = 500,
                                     qvalueCutoff = 1,
                                     pvalueCutoff = 1, pAdjustMethod = "BH")
  
  min(enrichKEGG_PL23vR848$pvalue)
  min(enrichKEGG_PL23vR848$p.adjust)
  nrow(subset(enrichKEGG_PL23vR848, enrichKEGG_PL23vR848$pvalue < 0.05))
  nrow(subset(enrichKEGG_PL23vR848, enrichKEGG_PL23vR848$p.adjust < 0.05))
  head(enrichKEGG_PL23vR848)
  subset(enrichKEGG_PL23vR848, enrichKEGG_PL23vR848$p.adjust < 0.05)$Description
  
  subset(enrichKEGG_PL23vR848, enrichKEGG_PL23vR848$Description == "JAK-STAT signaling pathway - Mus musculus (house mouse)")
  
# KEGG pathway gene set enrichment analysis ----------------------------------
  keggList_PL23 <- mk_list(subset(PL23_entrez, padj <= 0.05))
  gseKEGG_PL23 <- gseKEGG(geneList = keggList_PL23, 
                          organism = "mmu", keyType = "ncbi-geneid",
                          # minGSSize = 10, maxGSSize = 500, 
                          minGSSize = 5, maxGSSize = 500, 
                          nPerm = 10000,
                          pvalueCutoff = 1, pAdjustMethod = "BH")
  min(gseKEGG_PL23$pvalue)
  min(gseKEGG_PL23$p.adjust)
  nrow(subset(gseKEGG_PL23, gseKEGG_PL23$pvalue < 0.05))
  nrow(subset(gseKEGG_PL23, gseKEGG_PL23$p.adjust < 0.05))
  head(gseKEGG_PL23)
  gseKEGG_PL23[c(1:10), ]

  # Make dot plot --------------------------------------------------------------
    dotplot(gseKEGG_PL23, showCategory = 10, 
            color = "NES", x = "p.adjust", size = "Count") +
    ggtitle("PL2-3 + 2DG vs PL2-3\nKEGG GSEA (Top 10 Enriched Pathways)")
  
  keggList_PL23vR848 <- mk_list(subset(PL23vR848_entrez, padj <= 0.05))
  gseKEGG_PL23vR848 <- gseKEGG(geneList = keggList_PL23vR848, 
                          organism = "mmu", keyType = "ncbi-geneid",
                          minGSSize = 10, maxGSSize = 500, 
                          nPerm = 10000,
                          pvalueCutoff = 1, pAdjustMethod = "none")
  
  min(gseKEGG_PL23vR848$pvalue)
  min(gseKEGG_PL23vR848$p.adjust)
  nrow(subset(gseKEGG_PL23vR848, gseKEGG_PL23vR848$pvalue < 0.05))
  nrow(subset(gseKEGG_PL23vR848, gseKEGG_PL23vR848$p.adjust < 0.05))
  head(gseKEGG_PL23vR848)
  
  # Make dot plot --------------------------------------------------------------
    dotplot(gseKEGG_PL23vR848, showCategory = 10, 
            color = "NES", x = "p.adjust", size = "Count") +
    ggtitle("PL2-3 vs R848\nKEGG GSEA (Top 10 Enriched Pathways)")
  
  
# Pathview visualization -------------------------------------------------------
  # # JAK-STAT signaling pathway - Mus musculus (house mouse)
  #   pathview(gene.data = keggList_PL23vR848, 
  #            pathway.id = "mmu04630",
  #            species = "mmu", 
  #            gene.idtype = "entrez", 
  #            kegg.native = TRUE,
  #            expand.node = TRUE,
  #            split.group = TRUE,
  #            limit = list(gene = max(abs(keggList_PL23)), cpd = 1),
  #            # bins = list(gene = 10, cpd = 10), 
  #            # both.dirs = list(gene = T, cpd = T), 
  #            # low = list(gene = "green", cpd = "blue"), 
  #            # mid = list(gene = "gray", cpd = "gray"), 
  #            # high = list(gene = "red", cpd = "yellow"), 
  #            na.col = "transparent")  
  
  
  # # Metabolic pathways - Mus musculus (house mouse)
  # pathview(gene.data = keggList_PL23, gene.idtype = "entrez",
  #          pathway.id = 'mmu01100',
  #          species = "mmu",
  #          limit = list(gene=max(abs(keggList_PL23)), cpd=1))
  
  
  # Systemic lupus erythematosus - Mus musculus (house mouse)
  pathview(gene.data = keggList_PL23, gene.idtype = "entrez",
           pathway.id = 'mmu05322',
           species = "mmu",
           limit = list(gene=max(abs(keggList_PL23)), cpd=1))
  # 
  # # keggview.native(plot.data.gene = keggList_PL23, gene.idtype = "entrez",
  # #                 pathway.name = "mmu05322",
  # #                 species = "mmu",
  # #                 limit = list(gene=max(abs(keggList_PL23)), cpd=1))
  # 
  # # pathview(gene.data = keggList_PL23, 
  # #          pathway.id = "mmu05322",
  # #          species = "mmu", 
  # #          gene.idtype = "entrez", 
  # #          expand.node = TRUE,
  # #          split.group = TRUE,
  # #          limit = list(gene = max(abs(keggList_PL23)), cpd = 1),
  # #          # bins = list(gene = 10, cpd = 10), 
  # #          # low = list(gene = "green", cpd = "blue"), 
  # #          # mid = list(gene = "gray", cpd = "gray"), 
  # #          # high = list(gene = "red", cpd = "yellow"), 
  # #          na.col = "transparent")  
           
           
  # # Glycolysis pathway
  # pathview(gene.data = keggList_PL23, gene.idtype = "entrez",
  #          pathway.id = 'mmu00010',
  #          species = "mmu",
  #          limit = list(gene=max(abs(keggList_PL23)), cpd=1))
  

  
  
# Reactome pathway over-representation analysis --------------------------------
  reactomePA_PL23 <- enrichPathway(gene = names(enrichList_PL23), 
                                   organism = "mouse",
                                   minGSSize = 5, maxGSSize = 200,
                                   pvalueCutoff = 1, pAdjustMethod = "BH", 
                                   qvalueCutoff = 1, readable = TRUE)
  min(reactomePA_PL23$pvalue)
  min(reactomePA_PL23$p.adjust)
  nrow(subset(reactomePA_PL23, reactomePA_PL23$pvalue < 0.05))
  nrow(subset(reactomePA_PL23, reactomePA_PL23$p.adjust < 0.05))
  head(reactomePA_PL23)
  
  subset(reactomePA_PL23, reactomePA_PL23$p.adjust < 0.05)
  subset(reactomePA_PL23, reactomePA_PL23$pvalue < 0.05)$Description
  
# Reactome pathway gene set enrichment analysis --------------------------------
  reactomeGSA_PL23 <- gsePathway(geneList = entrezList_PL23, 
                                 organism = "mouse",
                                 minGSSize = 5, maxGSSize = 100,
                                 pvalueCutoff = 1, pAdjustMethod = "BH")
  min(reactomeGSA_PL23$pvalue)
  min(reactomeGSA_PL23$p.adjust)
  nrow(subset(reactomeGSA_PL23, reactomeGSA_PL23$pvalue < 0.05))
  nrow(subset(reactomeGSA_PL23, reactomeGSA_PL23$p.adjust < 0.05))
  head(reactomeGSA_PL23)
  
  # Make dot plot --------------------------------------------------------------
  dotplot(reactomeGSA_PL23, showCategory = 10, 
          color = "NES", x = "p.adjust", size = "Count") +
    ggtitle("PL2-3 + 2DG vs PL2-3\nReactome GSEA (Top 10 Enriched Pathways)")
  
  
  entrezList_PL23vR848 <- mk_list(subset(PL23vR848_entrez, padj <= 0.05))
  entrezList_PL23vR848 <- gsePathway(geneList = entrezList_PL23vR848, 
                                     organism = "mouse",
                                     minGSSize = 5, maxGSSize = 100,
                                     pvalueCutoff = 1, pAdjustMethod = "BH")
  min(entrezList_PL23vR848$pvalue)
  min(entrezList_PL23vR848$p.adjust)
  nrow(subset(entrezList_PL23vR848, entrezList_PL23vR848$pvalue < 0.05))
  nrow(subset(entrezList_PL23vR848, entrezList_PL23vR848$p.adjust < 0.05))
  head(entrezList_PL23vR848)
  
  # Make dot plot --------------------------------------------------------------
  dotplot(entrezList_PL23vR848, showCategory = 10, 
          color = "NES", x = "p.adjust", size = "Count") +
    ggtitle("PL2-3 vs R848\nReactome GSEA (Top 10 Enriched Pathways)")
  
  
  
# Reactome Visualization -------------------------------------------------------
  viewPathway("Complement cascade", 
              organism = "mouse", foldChange = entrezList_PL23)

  viewPathway("Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell", 
              organism = "mouse", foldChange = entrezList_PL23)

  viewPathway("Immune System", 
              organism = "mouse", foldChange = entrezList_PL23)

  viewPathway("Innate Immune System", 
              organism = "mouse", foldChange = entrezList_PL23)
  
  viewPathway("Metabolism of lipids", organism = "mouse", foldChange = entrezList_PL23)
  