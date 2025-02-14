gsea_VP <- function(geneset_file, vpdata, plot_file_name, main_title, sub_title, xlims, ylims){
  geneset <- read.delim(geneset_file)
  geneset <- geneset[, c(2:7)]
  colnames(vpdata) <- c("SYMBOL", "log2FoldChange", "padj")
  geneset <- left_join(geneset, vpdata)
  if(missing(xlims)) xlims <- c(-10, 10)
  if(missing(ylims)) ylims <- c(0, 15)
  EnhancedVolcano(geneset, lab = geneset$SYMBOL, pCutoff = 0.05, FCcutoff = 1,
                  x = 'log2FoldChange', y = 'padj', 
                  xlim = xlims, ylim = ylims,
                  title = main_title, subtitle = sub_title, 
                  legendLabels = c("NS", expression(Log[2] ~ FC > 1), 
                                   expression(p - value < ~ 0.05),
                                   expression(p - value < ~ 0.05 ~ and ~ log[2] ~ FC > ~ 1)),
                  drawConnectors = TRUE, min.segment.length = 1, 
                  max.overlaps = 12, labSize = 4)
  ggsave(stri_join(c("./../../../Hallmark_VolcanoPlots/", plot_file_name), collapse = ""),
         units = "px", width = 3000, height = 3000, dpi = 300)
}


# Hallmark Inflammatory Response -----------------------------------------------
  head(hallmark$PL23_2DG_v_Ctrl)
  
  HALLMARK_INFLAMMATORY_RESPONSE_PL23 <- HALLMARK_INFLAMMATORY_RESPONSE_PL23[, c(2:7)]
  nrow(HALLMARK_INFLAMMATORY_RESPONSE_PL23)
  HALLMARK_INFLAMMATORY_RESPONSE_PL23 <- left_join(HALLMARK_INFLAMMATORY_RESPONSE_PL23, vpData_PL23)
  nrow(HALLMARK_INFLAMMATORY_RESPONSE_PL23)
  head(HALLMARK_INFLAMMATORY_RESPONSE_PL23)
  
  
# PL2-3 + 2DG vs PL2-3 ---------------------------------------------------------
  setwd("./handling_GSEA_results/GSEA_results/Hallmark_PL23_2DG_vs_PL23/gene_sets")

  max(vpData_PL23$log2FoldChange)
  min(vpData_PL23$log2FoldChange)
  -1*log10(min(vpData_PL23$padj))
  
  gsea_VP("HALLMARK_INFLAMMATORY_RESPONSE.tsv", vpData_PL23, 
          "INFLAMMATORY_RESPONSE_PL23_2DG_vs_PL23.png", 
          "Hallmark: Inflammatory Response", "PL2-3 + 2DG vs PL2-3",
          xlims = c(-10, 10), ylims = c(0, 4))
  gsea_VP("HALLMARK_INTERFERON_GAMMA_RESPONSE.tsv", vpData_PL23, 
          "INTERFERON_GAMMA_RESPONSE_PL23_2DG_vs_PL23.png", 
          "Hallmark: INTERFERON_GAMMA_RESPONSE", "PL2-3 + 2DG vs PL2-3",
          xlims = c(-10, 6), ylims = c(0, 7))
  gsea_VP("HALLMARK_IL2_STAT5_SIGNALING.tsv", vpData_PL23, 
          "IL2_STAT5_SIGNALING_PL23_2DG_vs_PL23.png", 
          "Hallmark: IL2_STAT5_SIGNALING", "PL2-3 + 2DG vs PL2-3",
          xlims = c(-10, 6), ylims = c(0, 7))
  gsea_VP("HALLMARK_IL6_JAK_STAT3_SIGNALING.tsv", vpData_PL23, 
          "IL6_JAK_STAT3_SIGNALING_PL23_2DG_vs_PL23.png", 
          "Hallmark: IL6_JAK_STAT3_SIGNALING", "PL2-3 + 2DG vs PL2-3",
          xlims = c(-10, 6), ylims = c(0, 7))
  gsea_VP("HALLMARK_COMPLEMENT.tsv", vpData_PL23, 
          "COMPLEMENT_PL23_2DG_vs_PL23.png", 
          "Hallmark: COMPLEMENT", "PL2-3 + 2DG vs PL2-3",
          xlims = c(-10, 6), ylims = c(0, 7))
  gsea_VP("HALLMARK_GLYCOLYSIS.tsv", vpData_PL23, 
          "GLYCOLYSIS_PL23_2DG_vs_PL23.png", 
          "Hallmark: GLYCOLYSIS", "PL2-3 + 2DG vs PL2-3",
          xlims = c(-10, 6), ylims = c(0, 7))
  gsea_VP("HALLMARK_INTERFERON_ALPHA_RESPONSE.tsv", vpData_PL23, 
          "INTERFERON_ALPHA_RESPONSE_PL23_2DG_vs_PL23.png", 
          "Hallmark: INTERFERON_ALPHA_RESPONSE", "PL2-3 + 2DG vs PL2-3",
          xlims = c(-10, 6), ylims = c(0, 7))
  gsea_VP("HALLMARK_HYPOXIA.tsv", vpData_PL23, 
          "HYPOXIA_PL23_2DG_vs_PL23.png", 
          "Hallmark: HYPOXIA", "PL2-3 + 2DG vs PL2-3",
          xlims = c(-10, 6), ylims = c(0, 7))
  gsea_VP("HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY.tsv", vpData_PL23, 
          "REACTIVE_OXYGEN_SPECIES_PATHWAY_PL23_2DG_vs_PL23.png", 
          "Hallmark: REACTIVE_OXYGEN_SPECIES_PATHWAY", "PL2-3 + 2DG vs PL2-3",
          xlims = c(-10, 6), ylims = c(0, 7))
  gsea_VP("HALLMARK_FATTY_ACID_METABOLISM.tsv", vpData_PL23, 
          "FATTY_ACID_METABOLISM_PL23_2DG_vs_PL23.png", 
          "Hallmark: FATTY_ACID_METABOLISM", "PL2-3 + 2DG vs PL2-3",
          xlims = c(-10, 6), ylims = c(0, 7))
  gsea_VP("HALLMARK_OXIDATIVE_PHOSPHORYLATION.tsv", vpData_PL23, 
          "OXIDATIVE_PHOSPHORYLATION_PL23_2DG_vs_PL23.png", 
          "Hallmark: OXIDATIVE_PHOSPHORYLATION", "PL2-3 + 2DG vs PL2-3",
          xlims = c(-10, 6), ylims = c(0, 7))
  gsea_VP("HALLMARK_G2M_CHECKPOINT.tsv", vpData_PL23, 
          "G2M_CHECKPOINT_PL23_2DG_vs_PL23.png", 
          "Hallmark: G2M_CHECKPOINT", "PL2-3 + 2DG vs PL2-3",
          xlims = c(-10, 6), ylims = c(0, 7))
  gsea_VP("HALLMARK_MYC_TARGETS_V1.tsv", vpData_PL23, 
          "MYC_TARGETS_V1_PL23_2DG_vs_PL23.png", 
          "Hallmark: MYC_TARGETS_V1", "PL2-3 + 2DG vs PL2-3",
          xlims = c(-10, 6), ylims = c(0, 7))
  gsea_VP("HALLMARK_E2F_TARGETS.tsv", vpData_PL23, 
          "E2F_TARGETS_PL23_2DG_vs_PL23.png", 
          "Hallmark: E2F_TARGETS", "PL2-3 + 2DG vs PL2-3",
          xlims = c(-10, 6), ylims = c(0, 7))
  

  
  
  
  
  
  
  # reset working directory ----------------------------------------------------
    setwd("./../../../..")
    getwd()
    
# R848 + 2DG vs R848 ---------------------------------------------------------
  setwd("./handling_GSEA_results/GSEA_results/Hallmark_R848_2DG_vs_R848/gene_sets")
  
  HALLMARK_INFLAMMATORY_RESPONSE_R848 <- read.delim("HALLMARK_INFLAMMATORY_RESPONSE.tsv")
  head(HALLMARK_INFLAMMATORY_RESPONSE_R848)
  
  gsea_VP("HALLMARK_INFLAMMATORY_RESPONSE.tsv", vpData_R848, 
          "INFLAMMATORY_RESPONSE_R848_2DG_vs_R848.png", 
          "Hallmark: Inflammatory Response", "R848 + 2DG vs R848",
          xlims = c(-10, 10), ylims = c(0, 4))
  
  # reset working directory ----------------------------------------------------
    setwd("./../../../..")
    getwd()
    