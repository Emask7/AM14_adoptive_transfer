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

# PL2-3: 2DG vs Control --------------------------------------------------------
  # Set up DESeq Data Set ------------------------------------------------------
    dds_PL23 <- DESeqDataSetFromMatrix(countData = cts[, 1:10], 
                                       colData = coldata[1:10, ],
                                       design = ~Drug)
    dds_PL23
    
    dds_PL23$Drug <- relevel(dds_PL23$Drug, ref = "Control")
    keep <- rowSums(counts(dds_PL23) >= 10) >= 5
    dds_PL23 <- dds_PL23[keep,]
    dds_PL23 <- DESeq(dds_PL23)
    resultsNames(dds_PL23)
    
  # Transform data -------------------------------------------------------------
    vsd <- vst(dds_PL23)
    rld <- rlog(dds_PL23)
    ntd <- normTransform(dds_PL23)
    
  # Heatmap of count matrix ----------------------------------------------------
    select <- order(rowMeans(counts(dds_PL23,normalized=TRUE)), 
                    decreasing=TRUE)[1:20]
    df <- as.data.frame(colData(dds_PL23)[, c("Treatment", "Cohort")])
    
    ntd_heatmap <- pheatmap(assay(ntd)[select,], cluster_rows=FALSE,
                            show_rownames=FALSE, cluster_cols=TRUE,
                            annotation_col=df,
                            labels_col = colData(dds_PL23)$Label_Name,
                            main = "Normalized Counts Transformation")
    
    vsd_heatmap <- pheatmap(assay(vsd)[select,], cluster_rows=FALSE,
                            show_rownames=FALSE, cluster_cols=TRUE,
                            annotation_col=df,
                            labels_col = colData(dds_PL23)$Label_Name,
                            main = "Variance Stabilizing Transformation")
    
    rld_heatmap <- pheatmap(assay(rld)[select,], cluster_rows=FALSE,
                            show_rownames=FALSE, cluster_cols=TRUE,
                            annotation_col=df,
                            labels_col = colData(dds_PL23)$Label_Name,
                            main = "Regularized Log Transformation")
    
  # Heatmap of sample-to-sample distances --------------------------------------
    sampleDists <- dist(t(assay(vsd)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(vsd$Label_Name)
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    sampleDist_heatmap <- pheatmap(sampleDistMatrix,
                                   clustering_distance_rows=sampleDists,
                                   clustering_distance_cols=sampleDists,
                                   col=colors)
    
  # Principal component plot ---------------------------------------------------
    pcaData <- plotPCA(vsd, intgroup=c("Treatment", "Cohort"), returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    ggplot(pcaData, aes(PC1, PC2, color=Treatment, shape=Cohort)) +
      geom_point(size=3) +
      geom_text(label=coldata$Mouse_Number[1:10], nudge_x=1, nudge_y=0, check_overlap=F) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      coord_fixed() +
      labs(title = "PL2-3 + 2DG vs PL2-3\nBefore correcting for batch effects")
    
  # PCA plot removing batch effects --------------------------------------------
    mat <- assay(vsd)
    mm <- model.matrix(~Stim + Drug + Stim:Drug, colData(vsd))
    mat <- removeBatchEffect(mat, batch=vsd$Cohort, design=mm)
    assay(vsd) <- mat
    pcaData_afterBatch <- plotPCA(vsd, intgroup=c("Treatment", "Cohort"), returnData=TRUE)
    percentVar_afterBatch <- round(100 * attr(pcaData, "percentVar"))
    ggplot(pcaData_afterBatch, aes(PC1, PC2, color=Treatment, shape=Cohort)) +
      geom_point(size=3) +
      geom_text(label=coldata$Mouse_Number[1:10], nudge_x=1, nudge_y=0, check_overlap=F) +
      xlab(paste0("PC1: ",percentVar_afterBatch[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar_afterBatch[2],"% variance")) +
      coord_fixed() +
      labs(title = "PL2-3 + 2DG vs PL2-3\nAfter correcting for batch effects")
    
  # Differential Expression Analysis -------------------------------------------
    res_PL23 <- results(dds_PL23)
    res_PL23
    summary(res_PL23, alpha = 0.05)
    nrow(subset(res_PL23, 
                res_PL23$padj <= 0.05 & res_PL23$log2FoldChange >= 1))
    nrow(subset(res_PL23, 
                res_PL23$padj <= 0.05 & res_PL23$log2FoldChange <= -1))
    
  # LFC shrinkage results ----------------------------------------------------
    res_lfcShrink_PL23 <- lfcShrink(dds_PL23, 
                                    coef = "Drug_2DG_vs_Control", 
                                    type = "apeglm")
    res_lfcShrink_PL23
    summary(res_lfcShrink_PL23, alpha = 0.05)
    nrow(subset(res_lfcShrink_PL23, 
                res_lfcShrink_PL23$padj <= 0.05 & 
                  res_lfcShrink_PL23$log2FoldChange >= 1))
    nrow(subset(res_lfcShrink_PL23, 
                res_lfcShrink_PL23$padj <= 0.05 & 
                  res_lfcShrink_PL23$log2FoldChange <= -1))
