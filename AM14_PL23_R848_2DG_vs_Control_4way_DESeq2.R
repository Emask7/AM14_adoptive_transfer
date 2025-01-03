# Import count data and sample info --------------------------------------------
  cts <- read.delim("raw_data/AM14_rawCounts_No_Ighv_or_Igkv_genes.txt")
  colnames(cts)
  cts <- cts[, c(2, 5:24)]
  head(cts)
  nrow(cts)
  cts <- cts[!duplicated(cts$external_gene_name), ]
  nrow(cts)
  rownames(cts) <- cts$external_gene_name
  cts <- subset(cts, external_gene_name != "7SK" & 
                  external_gene_name != "5_8S_rRNA" & 
                  external_gene_name != "5S_rRNA")
  nrow(cts)
  cts <- as.matrix(cts[, 2:21])
  head(cts)

  coldata <- read.csv("sample_info.csv")
  coldata <- coldata[1:20 , c(1:3, 5:6)]
  coldata$Cohort <- factor(coldata$Cohort)
  coldata$Treatment <- factor(coldata$Treatment)
  coldata
  
# Make sure the columns of cts and rows of coldata are in the same order -------
  for (x in 1:ncol(cts)) {
    if(colnames(cts)[x] == coldata$Sample[x]) print(c(x, "true")) 
    else print(c(x, "false"))
  }

# Set up DESeq Data Set --------------------------------------------------------
  dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata,
                                design = ~Cohort + Treatment)
  dds
  keep <- rowSums(counts(dds) >= 10) >= 5
  dds <- dds[keep,]
  dds <- DESeq(dds)
  resultsNames(dds)
  
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
  
# Differential Expression Analysis -------------------------------------------
  resultsNames(dds)
    
  res_PL23 <- results(dds, contrast = c("Treatment", "PL2-3+2DG", "PL2-3"))
  res_PL23 <- results(dds, contrast = c("Treatment", "PL2-3+2DG", "PL2-3"), 
                      lfcThreshold = 0.6, alpha = 0.05)
  res_PL23
  summary(res_PL23)
  nrow(subset(res_PL23, 
              res_PL23$padj <= 0.05 & res_PL23$log2FoldChange >= 0.6))
  nrow(subset(res_PL23, 
              res_PL23$padj <= 0.05 & res_PL23$log2FoldChange <= -0.6))
  
  rownames(subset(res_PL23, res_PL23$padj <= 0.05 & res_PL23$log2FoldChange <= -0.6))
  
  res_R848 <- results(dds, contrast = c("Treatment", "R848+2DG", "R848"))
  res_R848 <- results(dds, contrast = c("Treatment", "R848+2DG", "R848"), 
                      lfcThreshold = 0.6, alpha = 0.05)
  res_R848
  summary(res_R848)
  nrow(subset(res_R848, 
              res_R848$padj <= 0.05 & res_R848$log2FoldChange >= 1))
  nrow(subset(res_R848, 
              res_R848$padj <= 0.05 & res_R848$log2FoldChange <= -1))
  
  res_PL23vR848 <- results(dds, contrast = c("Treatment", "PL2-3", "R848"))
  res_PL23vR848 <- results(dds, contrast = c("Treatment", "PL2-3", "R848"), 
                           lfcThreshold = 0.6, alpha = 0.05)
  res_PL23vR848
  summary(res_PL23vR848)
  nrow(subset(res_PL23vR848, 
              res_PL23vR848$padj <= 0.05 & res_PL23vR848$log2FoldChange >= 1))
  nrow(subset(res_PL23vR848, 
              res_PL23vR848$padj <= 0.05 & res_PL23vR848$log2FoldChange <= -1))
  
  res_2DG_PL23vR848 <- results(dds, contrast = c("Treatment", "PL2-3+2DG", "R848+2DG"))
  res_2DG_PL23vR848 <- results(dds, contrast = c("Treatment", "PL2-3+2DG", "R848+2DG"), 
                               lfcThreshold = 1, alpha = 0.05)
  res_2DG_PL23vR848
  summary(res_2DG_PL23vR848)
  nrow(subset(res_2DG_PL23vR848, 
              res_2DG_PL23vR848$padj <= 0.05 & res_2DG_PL23vR848$log2FoldChange >= 1))
  nrow(subset(res_2DG_PL23vR848, 
              res_2DG_PL23vR848$padj <= 0.05 & res_2DG_PL23vR848$log2FoldChange <= -1))
  
  pl23_down <- rownames(data.frame(subset(res_PL23, 
                       res_PL23$padj <= 0.05 & res_PL23$log2FoldChange <= -1)))
  pl23_down <- data.frame(pl23_down, pl23_down)
  colnames(pl23_down) <- c("gene", "pl23_2dg_ctrl")
  
  pl23_r848_down <- rownames(data.frame(subset(res_PL23vR848, 
                           res_PL23vR848$padj <= 0.05 & res_PL23vR848$log2FoldChange <= -1)))
  pl23_r848_down <- data.frame(pl23_r848_down, pl23_r848_down)
  colnames(pl23_r848_down) <- c("gene", "pl23_r848")
  
  temp_join <- inner_join(pl23_down, pl23_r848_down, by = "gene")
  head(temp_join)
  nrow(temp_join)
  
  # Filter out rows with a padj value of NA ------------------------------------
    # Note: There are a few reasons why a p value or padj value would be NA
    # According to the DESeq2 manual, these are the reasons:
      # If within a row, all samples have zero counts, the baseMean column will
        # be zero, and the LFC estimates, p value and padj will all be NA.
      # If a row contains a sample with an extreme count outlier then the 
        # p value and padj will be set to NA.
      # If a row is filtered by automatic independent filtering, for having a 
        # low mean normalized count, then only padj will be set to NA.
  
    PL23_filt_res <- data.frame(subset(res_PL23, !is.na(padj)))
    R848_filt_res <- data.frame(subset(res_R848, !is.na(padj)))
    PL23vR848_filt_res <- data.frame(subset(res_PL23vR848, !is.na(padj)))

# Save count data to Excel file ------------------------------------------------
  wb <- createWorkbook("Output/DESeq2_counts.xlsx")

  addWorksheet(wb, "PL23_2DG_vs_Ctrl")
  addWorksheet(wb, "R848_2DG_vs_Ctrl")
  addWorksheet(wb, "PL23_Ctrl_vs_R848_Ctrl")

  writeData(wb, "PL23_2DG_vs_Ctrl", PL23_filt_res, rowNames = TRUE)
  writeData(wb, "R848_2DG_vs_Ctrl", R848_filt_res, rowNames = TRUE)
  writeData(wb, "PL23_Ctrl_vs_R848_Ctrl", PL23vR848_filt_res, rowNames = TRUE)
  
  saveWorkbook(wb, "Output/DESeq2_counts.xlsx", overwrite = TRUE)