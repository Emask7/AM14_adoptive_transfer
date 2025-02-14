# Functions --------------------------------------------------------------------
  summary_v2 <- function(res, title, p_cutoff, lfc_cutoff){
    if(missing(p_cutoff)) p_cutoff <- 0.05
    if(missing(lfc_cutoff)) lfc_cutoff <- 0.6
    if(missing(title)) title <- "# Genes"
    
    up_string <- stri_join(c("Up (LFC >= ", lfc_cutoff, ")"), collapse = "")
    down_string <- stri_join(c("Down (LFC <= -", lfc_cutoff, ")"), collapse = "")
    
    DEG_summary <- data.frame(c(up_string, down_string),
                              c(nrow(subset(res, res$padj <= p_cutoff & res$log2FoldChange >= lfc_cutoff)),
                                nrow(subset(res, res$padj <= p_cutoff & res$log2FoldChange <= (-1*lfc_cutoff)))))
    
    col1_name <- stri_join(c("padj <= ", p_cutoff), collapse = "")
    colnames(DEG_summary) <- c(col1_name, title)
    DEG_summary
  }

  write_DEG_CSV <- function(res, lfc_cutoff, file_start){
    up <- subset(res, res$padj <= 0.05 & res$log2FoldChange >= lfc_cutoff)[, c(1, 2, 4)]
    down <- subset(res, res$padj <= 0.05 & res$log2FoldChange <= (-1*lfc_cutoff))[, c(1, 2, 4)]
    all <- subset(res, res$padj <= 0.05 & abs(res$log2FoldChange) >= lfc_cutoff)[, c(1, 2, 4)]
    
    if(!is.null(file_start)){
      write.table(up, 
                  file = stri_join(c("Output/Gene_Lists/", file_start, "_Up.csv"), collapse = ""),
                  sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
      write.table(down, 
                  file = stri_join(c("Output/Gene_Lists/", file_start, "_Down.csv"), collapse = ""),
                  sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
      write.table(all,
                  file = stri_join(c("Output/Gene_Lists/", file_start, "_All.csv"), collapse = ""),
                  sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
    }
    
    return(list(up = up, down = down, all = all))
  }

  write_sig_LFCs <- function(res, lfc_cutoff, ID_type, file_start){
    dat <- subset(res, res$padj <= 0.05 & abs(res$log2FoldChange) >= lfc_cutoff)
    
    if(missing(ID_type)) print("specify ID type: ensembl_gene_id, external_gene_name, or entrezgene")
    else if(ID_type == "ensembl_gene_id") dat <- dat[, c(1, 6)]
    else if(ID_type == "external_gene_name") dat <- dat[, c(2, 6)]
    else if(ID_type == "entrezgene") dat <- dat[, c(4, 6)]
    else print("specify ID type: ensembl_gene_id, external_gene_name, or entrezgene")
    
    colnames(dat) <- c("#", "LogFoldChange")
    
    if(!is.null(file_start)){
      write.table(dat,
                  file = stri_join(c("Output/Gene_Lists/", file_start, "_LFC_values.csv"), collapse = ""),
                  sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
    }
    
    dat
  }
  

# Import count data and sample info --------------------------------------------
  cts <- read.delim("raw_data/AM14_rawCounts_No_Ighv_or_Igkv_genes.txt")
  colnames(cts)
  nrow(cts)
  cts <- subset(cts, external_gene_name != "7SK" & 
                  external_gene_name != "5_8S_rRNA" & 
                  external_gene_name != "5S_rRNA")
  nrow(cts)
  
  # # For gene symbols as row names ----------------------------------------------
  #   cts <- cts[!duplicated(cts$external_gene_name), ]
  #   nrow(cts)
  #   cts <- cts[!is.na(cts$external_gene_name), ]
  #   nrow(cts)
  #   rownames(cts) <- cts$external_gene_name
  #   nrow(cts)

  # For Ensembl IDs as row names -----------------------------------------------
    cts <- cts[!duplicated(cts$ensembl_gene_id), ]
    nrow(cts)
    cts <- cts[!is.na(cts$ensembl_gene_id), ]
    nrow(cts)
    rownames(cts) <- cts$ensembl_gene_id
    nrow(cts)
    
    gene_IDs <- cts[, c(1:4)]
    head(gene_IDs)

  # # For Entrez IDs as row names ------------------------------------------------
  #   cts <- cts[!duplicated(cts$entrezgene), ]
  #   nrow(cts)
  #   cts <- cts[!is.na(cts$entrezgene), ]
  #   nrow(cts)
  #   rownames(cts) <- cts$entrezgene
  #   nrow(cts)

  # Remove unneeded columns and set up experimental factors --------------------
    cts <- as.matrix(cts[, 5:24])
    head(cts)
  
    coldata <- read.csv("sample_info.csv")
    coldata <- coldata[1:20 , c(1:3, 5:6)]
    coldata$Cohort <- factor(coldata$Cohort)
    coldata$Treatment <- factor(coldata$Treatment)
    coldata
  
    # Make sure the columns of cts and rows of coldata are in the same order ---
      for (x in 1:ncol(cts)) {
        if(colnames(cts)[x] == coldata$Sample[x]) print(c(x, "true")) 
        else print(c(x, "false"))
      }
      rm(x)

# Set up DESeq Data Set --------------------------------------------------------
  dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata,
                                design = ~Cohort + Treatment)
  dds
  keep <- rowSums(counts(dds) >= 10) >= 5
  dds <- dds[keep,]
  dds <- DESeq(dds)
  resultsNames(dds)
  rm(keep)
  
# Quality Check Steps ----------------------------------------------------------
  # Transform data -------------------------------------------------------------
    vsd <- vst(dds)
    rld <- rlog(dds)
    ntd <- normTransform(dds)

  # Heatmap of count matrix ----------------------------------------------------
    select <- order(rowMeans(counts(dds,normalized=TRUE)), 
                    decreasing=TRUE)[1:20]
    df <- as.data.frame(colData(dds)[, c("Treatment", "Cohort")])
    
    png(filename = "Output/Heatmaps/Heatmap - Normalized Counts Transformation.png", 
        width = 1000, height = 1000, units = "px", pointsize = 10, res = 200, 
        bg = "white", family = "", type = "windows", symbolfamily="default")
    pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE, 
             cluster_cols=TRUE, annotation_col=df,
             labels_col = colData(dds)$Label_Name,
             main = "Normalized Counts Transformation")
    dev.off()
    
    png(filename = "Output/Heatmaps/Heatmap - Variance Stabilizing Transformation.png", 
        width = 1000, height = 1000, units = "px", pointsize = 10, res = 200, 
        bg = "white", family = "", type = "windows", symbolfamily="default")
    pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE, 
             cluster_cols=TRUE, annotation_col=df,
             labels_col = colData(dds)$Label_Name,
             main = "Variance Stabilizing Transformation")
    dev.off()
    
    png(filename = "Output/Heatmaps/Heatmap - Regularized Log Transformation.png", 
        width = 1000, height = 1000, units = "px", pointsize = 10, res = 200, 
        bg = "white", family = "", type = "windows", symbolfamily="default")
    pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE, 
             cluster_cols=TRUE, annotation_col=df,
             labels_col = colData(dds)$Label_Name,
             main = "Regularized Log Transformation")
    dev.off()
    
  # Heatmap of sample-to-sample distances --------------------------------------
    sampleDists <- dist(t(assay(vsd)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(vsd$Label_Name)
    colnames(sampleDistMatrix) <- NULL
  
    png(filename = "Output/Heatmaps/Heatmap - Sample-to-Sample Distances.png", 
        width = 1200, height = 1200, units = "px", pointsize = 10, res = 200, 
        bg = "white", family = "", type = "windows", symbolfamily="default")
    pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists,
             clustering_distance_cols=sampleDists,
             col=colorRampPalette( rev(brewer.pal(9, "Blues")) )(255),
             main = "Sample-to-Sample Distances")
    dev.off()
    
  # Principal component plot ---------------------------------------------------
    pcaData <- plotPCA(vsd, intgroup=c("Treatment", "Cohort"), returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    
    png(filename = "Output/PCA Plot - Before Batch Correction.png", 
        width = 1500, height = 1500, units = "px", pointsize = 10, res = 200, 
        bg = "white", family = "", type = "windows", symbolfamily="default")
    ggplot(pcaData, aes(PC1, PC2, color=Treatment, shape=Cohort)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      coord_fixed() +
      labs(title = "Before correcting for batch effects")
    dev.off()
    
  # PCA plot removing batch effects --------------------------------------------
    mat <- assay(vsd)
    mm <- model.matrix(~Treatment, colData(vsd))
    mat <- removeBatchEffect(mat, batch=vsd$Cohort, design=mm)
    assay(vsd) <- mat
    pcaData <- plotPCA(vsd, intgroup=c("Treatment", "Cohort"), returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    
    png(filename = "Output/PCA Plot - After Batch Correction.png", 
        width = 1500, height = 1500, units = "px", pointsize = 10, res = 200, 
        bg = "white", family = "", type = "windows", symbolfamily="default")
    ggplot(pcaData, aes(PC1, PC2, color=Treatment, shape=Cohort)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      coord_fixed() +
      labs(title = "After correcting for batch effects")
    dev.off()
  
# Differential Expression Analysis ---------------------------------------------
  resultsNames(dds)
    
# Get results (default methods) ------------------------------------------------
  res_PL23 <- results(dds, contrast = c("Treatment", "PL2-3+2DG", "PL2-3"))
  res_PL23 <- data.frame(subset(res_PL23, !is.na(padj)))
  summary_v2(res_PL23, "PL2-3+2DG vs PL2-3")
  
  res_R848 <- results(dds, contrast = c("Treatment", "R848+2DG", "R848"))
  res_R848 <- data.frame(subset(res_R848, !is.na(padj)))
  summary_v2(res_R848, "R848+2DG vs R848")
  
  res_PL23vR848 <- results(dds, contrast = c("Treatment", "PL2-3", "R848"))
  res_PL23vR848 <- data.frame(subset(res_PL23vR848, !is.na(padj)))
  summary_v2(res_PL23vR848, "PL2-3 vs R848")
  
  res_2DG_PL23vR848 <- results(dds, contrast = c("Treatment", "PL2-3+2DG", "R848+2DG"))
  res_2DG_PL23vR848 <- data.frame(subset(res_2DG_PL23vR848, !is.na(padj)))
  summary_v2(res_2DG_PL23vR848, "PL2-3+2DG vs R848+2DG")
  
    # Note: There are a few reasons why a p value or padj value would be NA
    # According to the DESeq2 manual, these are the reasons:
    # If within a row, all samples have zero counts, the baseMean column will
    # be zero, and the LFC estimates, p value and padj will all be NA.
    # If a row contains a sample with an extreme count outlier then the 
    # p value and padj will be set to NA.
    # If a row is filtered by automatic independent filtering, for having a 
    # low mean normalized count, then only padj will be set to NA.
    
    
  # Make DEG lists and export to CSV files -------------------------------------
    write_DEG_CSV(res_PL23, 1, "PL2-3_2DG_vs_Control")
    write_DEG_CSV(res_R848, 1, "R848_2DG_vs_Control")
    write_DEG_CSV(res_PL23vR848, 1, "PL2-3_Ctrl_vs_R848_Ctrl")
    write_DEG_CSV(res_2DG_PL23vR848, 1, "PL2-3_2DG_vs_R848_2DG")

# # Get results (Threshold-Based Wald tests) -----------------------------------
#   res_PL23_lfc <- results(dds, contrast = c("Treatment", "PL2-3+2DG", "PL2-3"),
#                           lfcThreshold = 0.6, altHypothesis = "greaterAbs", alpha = 0.05)
#   res_R848_lfc <- results(dds, contrast = c("Treatment", "R848+2DG", "R848"),
#                           lfcThreshold = 0.6, altHypothesis = "greaterAbs", alpha = 0.05)
#   res_PL23vR848_lfc <- results(dds, contrast = c("Treatment", "PL2-3", "R848"),
#                                lfcThreshold = 0.6, altHypothesis = "greaterAbs", alpha = 0.05)
#   res_2DG_PL23vR848_lfc <- results(dds, contrast = c("Treatment", "PL2-3+2DG", "R848+2DG"),
#                                    lfcThreshold = 1, altHypothesis = "greaterAbs", alpha = 0.05)
# 
#   summary(res_PL23_lfc)
#   summary(res_R848_lfc)
#   summary(res_PL23vR848_lfc)
#   summary(res_2DG_PL23vR848_lfc)

# Save count data to Excel file ------------------------------------------------
  # Add columns with alternative gene identifiers to the DESeq2 results --------
    res_PL23_full <- tibble::rownames_to_column(res_PL23, var = "ensembl_gene_id")
    res_PL23_full <- right_join(gene_IDs, res_PL23_full)
    head(res_PL23_full)
    
    res_R848_full <- tibble::rownames_to_column(res_R848, var = "ensembl_gene_id")
    res_R848_full <- right_join(gene_IDs, res_R848_full)
    head(res_R848_full)
    
    res_PL23vR848_full <- tibble::rownames_to_column(res_PL23vR848, var = "ensembl_gene_id")
    res_PL23vR848_full <- right_join(gene_IDs, res_PL23vR848_full)
    head(res_PL23vR848_full)
    
    res_2DG_PL23vR848_full <- tibble::rownames_to_column(res_2DG_PL23vR848, var = "ensembl_gene_id")
    res_2DG_PL23vR848_full <- right_join(gene_IDs, res_2DG_PL23vR848_full)
    head(res_2DG_PL23vR848_full)
    
  # Make a data.frame that summarizes the numbers of DEGs detected -------------
    full_summary <- full_join(summary_v2(res_PL23, "PL2-3+2DG vs PL2-3"), 
                                summary_v2(res_PL23, "PL2-3+2DG vs PL2-3", lfc_cutoff = 1))
    
    full_summary_2 <- full_join(summary_v2(res_R848, "R848+2DG vs R848"), 
                                summary_v2(res_R848, "R848+2DG vs R848", lfc_cutoff = 1))
    
    full_summary_3 <- full_join(summary_v2(res_PL23vR848, "PL2-3 vs R848"), 
                                summary_v2(res_PL23vR848, "PL2-3 vs R848", lfc_cutoff = 1))
    
    full_summary_4 <- full_join(summary_v2(res_2DG_PL23vR848, "PL2-3+2DG vs R848+2DG"), 
                                summary_v2(res_2DG_PL23vR848, "PL2-3+2DG vs R848+2DG", lfc_cutoff = 1))
    
    full_summary <- full_join(full_summary, full_summary_2)
    full_summary <- full_join(full_summary, full_summary_3)
    full_summary <- full_join(full_summary, full_summary_4)
    full_summary
    
    rm(full_summary_2, full_summary_3, full_summary_4)
    
    
  # Save results to an Excel file ----------------------------------------------
    wb <- createWorkbook("Output/DESeq2_results.xlsx")
    
    addWorksheet(wb, "DEG Summaries")
    addWorksheet(wb, "PL23_2DG_vs_Ctrl")
    addWorksheet(wb, "R848_2DG_vs_Ctrl")
    addWorksheet(wb, "PL23_Ctrl_vs_R848_Ctrl")
    addWorksheet(wb, "PL23_2DG_vs_R848_2DG")
    addWorksheet(wb, "Raw_Gene_Counts")
    addWorksheet(wb, "Normalized_Gene_Counts")
    
    writeData(wb, "DEG Summaries", full_summary)
    writeData(wb, "PL23_2DG_vs_Ctrl", res_PL23_full)
    writeData(wb, "R848_2DG_vs_Ctrl", res_R848_full)
    writeData(wb, "PL23_Ctrl_vs_R848_Ctrl", res_PL23vR848_full)
    writeData(wb, "PL23_2DG_vs_R848_2DG", res_2DG_PL23vR848_full)
    writeData(wb, "Raw_Gene_Counts", 
              counts(dds, normalized = FALSE), rowNames = TRUE)
    writeData(wb, "Normalized_Gene_Counts", 
              counts(dds, normalized = TRUE), rowNames = TRUE)
    
    saveWorkbook(wb, "Output/DESeq2_results.xlsx", overwrite = TRUE)
    
# Make DEG Venn Diagram --------------------------------------------------------
  DEG_lists <- list(
    # PL23_2DGvsCtrl = write_DEG_CSV(res_PL23_full, 1, "PL2-3_2DG_vs_Control"),
    # R848_2DGvsCtrl = write_DEG_CSV(res_R848_full, 1, "R848_2DG_vs_Control"),
    # PL23vR848 = write_DEG_CSV(res_PL23vR848_full, 1, "PL2-3_Ctrl_vs_R848_Ctrl")
    PL23_2DGvsCtrl = write_sig_LFCs(res_PL23_full, 1, "external_gene_name", "PL2-3_2DG_vs_Control"),
    R848_2DGvsCtrl = write_sig_LFCs(res_R848_full, 1, "external_gene_name", "R848_2DG_vs_Control"),
    PL23vR848 = write_sig_LFCs(res_PL23vR848_full, 1, "external_gene_name", "PL2-3_Ctrl_vs_R848_Ctrl")
  )
  nrow(DEG_lists$PL23_2DGvsCtrl$down)
  nrow(DEG_lists$R848_2DGvsCtrl$down)
  nrow(DEG_lists$PL23vR848$down)
  
  head(DEG_lists$PL23_2DGvsCtrl)
  
  
  write_DEG_CSV(res_PL23, 1, "PL2-3_2DG_vs_Control")
  write_DEG_CSV(res_R848, 1, "R848_2DG_vs_Control")
  write_DEG_CSV(res_PL23vR848, 1, "PL2-3_Ctrl_vs_R848_Ctrl")
  write_DEG_CSV(res_2DG_PL23vR848, 1, "PL2-3_2DG_vs_R848_2DG")
  
    
  venn_up <- venndetail(list("PL2-3+2DG vs PL2-3" = DEG_lists$PL23_2DGvsCtrl$up$ensembl_gene_id,
                             "R848+2DG vs R848" = DEG_lists$R848_2DGvsCtrl$up$ensembl_gene_id,
                             "PL2-3 vs R848" = DEG_lists$PL23vR848$up$ensembl_gene_id))
  plot(venn_up, mycol = c("goldenrod1", "darkorange1", "red"), 
       filename = "Output/venn_diagram_upreg.png",
       margin = 0.1, cat.cex = 0.5)
  dev.off()
  detail(venn_up)
  

    
  venn_down <- venndetail(list("PL2-3+2DG vs PL2-3" = DEG_lists$PL23_2DGvsCtrl$down$ensembl_gene_id,
                               "R848+2DG vs R848" = DEG_lists$R848_2DGvsCtrl$down$ensembl_gene_id,
                               "PL2-3 vs R848" = DEG_lists$PL23vR848$down$ensembl_gene_id))
  plot(venn_down, mycol = c("darkseagreen1", "dodgerblue", "orchid"), 
       filename = "Output/venn_diagram_downreg.png",
       margin = 0.1, cat.cex = 0.5)
  dev.off()
  detail(venn_down)
