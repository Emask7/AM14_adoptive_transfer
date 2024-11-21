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
  coldata <- coldata[1:20 , c(1:3, 5:8)]
  coldata$Cohort <- factor(coldata$Cohort)
  coldata$Treatment <- factor(coldata$Treatment)
  coldata$Stim <- factor(coldata$Stim)
  coldata$Drug <- factor(coldata$Drug)
  coldata
  
# Make sure the columns of cts and rows of coldata are in the same order -------
  for (x in 1:ncol(cts)) {
    if(colnames(cts)[x] == coldata$Sample[x]) print(c(x, "true")) 
    else print(c(x, "false"))
  }

# All 4 Conditions -------------------------------------------------------------
  # Set up DESeq Data Set ------------------------------------------------------
    dds <- DESeqDataSetFromMatrix(countData = cts, 
                                  colData = coldata,
                                  design = ~Treatment + Cohort)
    dds
    keep <- rowSums(counts(dds) >= 10) >= 5
    dds <- dds[keep,]
    dds <- DESeq(dds)
    resultsNames(dds)
    
  # Transform data -----------------------------------------------------------
    vsd <- vst(dds)
    rld <- rlog(dds)
    ntd <- normTransform(dds)
    
  # Heatmap of count matrix --------------------------------------------------
    select <- order(rowMeans(counts(dds,normalized=TRUE)), 
                    decreasing=TRUE)[1:20]
    # df <- as.data.frame(colData(dds)[, c("Stim", "Drug", "Cohort")])
    df <- as.data.frame(colData(dds)[, c("Treatment", "Cohort")])
    
    ntd_heatmap <- pheatmap(assay(ntd)[select,], cluster_rows=FALSE,
                            show_rownames=FALSE, cluster_cols=TRUE,
                            annotation_col=df,
                            labels_col = colData(dds)$Label_Name,
                            main = "Normalized Counts Transformation")
    
    vsd_heatmap <- pheatmap(assay(vsd)[select,], cluster_rows=FALSE,
                            show_rownames=FALSE, cluster_cols=TRUE,
                            annotation_col=df,
                            labels_col = colData(dds)$Label_Name,
                            main = "Variance Stabilizing Transformation")
    
    rld_heatmap <- pheatmap(assay(rld)[select,], cluster_rows=FALSE,
                            show_rownames=FALSE, cluster_cols=TRUE,
                            annotation_col=df,
                            labels_col = colData(dds)$Label_Name,
                            main = "Regularized Log Transformation")
    
  # Heatmap of sample-to-sample distances ------------------------------------
    sampleDists <- dist(t(assay(vsd)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(vsd$Label_Name)
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    sampleDist_heatmap <- pheatmap(sampleDistMatrix,
                                   clustering_distance_rows=sampleDists,
                                   clustering_distance_cols=sampleDists,
                                   col=colors)
    
  # Principal component plot -------------------------------------------------
    pcaData <- plotPCA(vsd, intgroup=c("Treatment", "Cohort"), returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    ggplot(pcaData, aes(PC1, PC2, color=Treatment, shape=Cohort)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      coord_fixed() +
      labs(title = "Before correcting for batch effects")

  # PCA plot removing batch effects ------------------------------------------
    mat <- assay(vsd)
    mm <- model.matrix(~ ~Stim + Drug + Stim:Drug, colData(vsd))
    mat <- removeBatchEffect(mat, batch=vsd$Cohort, design=mm)
    assay(vsd) <- mat
    pcaData_afterBatch <- plotPCA(vsd, intgroup=c("Treatment", "Cohort"), returnData=TRUE)
    percentVar_afterBatch <- round(100 * attr(pcaData, "percentVar"))
    ggplot(pcaData_afterBatch, aes(PC1, PC2, color=Treatment, shape=Cohort)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar_afterBatch[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar_afterBatch[2],"% variance")) +
      coord_fixed() +
      labs(title = "After correcting for batch effects")
      
  # Differential Expression Analysis -------------------------------------------
    resultsNames(dds)
      
    res_PL23 <- results(dds, contrast = c("Treatment", "PL2-3+2DG", "PL2-3"))
    res_PL23
    summary(res_PL23, alpha = 0.05)
    nrow(subset(res_PL23, 
                res_PL23$padj <= 0.05 & res_PL23$log2FoldChange >= 1))
    nrow(subset(res_PL23, 
                res_PL23$padj <= 0.05 & res_PL23$log2FoldChange <= -1))
    
    res_R848 <- results(dds, contrast = c("Treatment", "R848+2DG", "R848"))
    res_R848
    summary(res_R848, alpha = 0.05)
    nrow(subset(res_R848, 
                res_R848$padj <= 0.05 & res_R848$log2FoldChange >= 1))
    nrow(subset(res_R848, 
                res_R848$padj <= 0.05 & res_R848$log2FoldChange <= -1))
    
    res_PL23vR848 <- results(dds, contrast = c("Treatment", "PL2-3", "R848"))
    res_PL23vR848
    summary(res_PL23vR848, alpha = 0.05)
    nrow(subset(res_PL23vR848, 
                res_PL23vR848$padj <= 0.05 & res_PL23vR848$log2FoldChange >= 1))
    nrow(subset(res_PL23vR848, 
                res_PL23vR848$padj <= 0.05 & res_PL23vR848$log2FoldChange <= -1))
      
  
      

# # Save count data to Excel file ------------------------------------------------
#   wb <- createWorkbook("DESeq2_counts.xlsx")
#   
#   addWorksheet(wb, "PL2-3_2DG_vs_Ctrl")
#   addWorksheet(wb, "R848_2DG_vs_Ctrl")
#   addWorksheet(wb, "PL23_all_vs_R848_all")
# 
#   writeData(wb, "PL2-3_2DG_vs_Ctrl", PL23.2DG_vs_PL23$dds)
#   
  