# Function that takes TSV results from all 4 comparisons and compiles them into 1 Excel file
  combine_results <- function(up, down){
    up <- read.delim(up)
    down <- read.delim(down)
    combined <- full_join(up, down)
    TLS <- separate_wider_delim(combined, cols = LEADING.EDGE, 
                                delim = stringr::regex("%|="),
                                names = c(NA, "TAGS", NA, "LIST", NA, "SIGNAL", NA))
    TLS <- data.frame(TLS$NAME, as.integer(TLS$TAGS), as.integer(TLS$LIST), as.integer(TLS$SIGNAL))
    colnames(TLS) <- c("FULL_NAME", "TAGS", "LIST", "SIGNAL")
    combined <- combined[, c(1:2, 4:10)]
    colnames(combined) <- c("FULL_NAME", "NAME", "SIZE", "ES", "NES", 
                            "NOM.p.val", "FDR.q.val", "FWER.p.val",  "RANK.AT.MAX")
    combined <- full_join(combined, TLS, by = "FULL_NAME")
    combined <- separate_wider_delim(combined, cols = NAME, 
                                     delim = stringr::regex("_"),
                                     names = c("SOURCE", "NAME"), too_many = "merge")
    combined <- combined[order(combined$FDR.q.val), ]
    return(data.frame(combined))
  }

  compile_results <- function(filename,
                              PL23_2DG_v_Ctrl_up, PL23_2DG_v_Ctrl_down,
                              R848_2DG_v_Ctrl_up, R848_2DG_v_Ctrl_down,
                              PL23_v_R848_up, PL23_v_R848_down,
                              PL23_2DG_v_R848_2DG_up, PL23_2DG_v_R848_2DG_down){
    
    PL23_2DG_v_Ctrl <- combine_results(PL23_2DG_v_Ctrl_up, PL23_2DG_v_Ctrl_down)
    R848_2DG_v_Ctrl <- combine_results(R848_2DG_v_Ctrl_up, R848_2DG_v_Ctrl_down)
    PL23_v_R848 <- combine_results(PL23_v_R848_up, PL23_v_R848_down)
    PL23_2DG_v_R848_2DG <- combine_results(PL23_2DG_v_R848_2DG_up, PL23_2DG_v_R848_2DG_down)
  
    if(!is.null(filename)){
      wb <- createWorkbook(filename)
    
      addWorksheet(wb, "PL23_2DG_vs_Ctrl")
      addWorksheet(wb, "R848_2DG_vs_Ctrl")
      addWorksheet(wb, "PL23_Ctrl_vs_R848_Ctrl")
      addWorksheet(wb, "PL23_2DG_vs_R848_2DG")
    
      writeData(wb, "PL23_2DG_vs_Ctrl", PL23_2DG_v_Ctrl)
      writeData(wb, "R848_2DG_vs_Ctrl", R848_2DG_v_Ctrl)
      writeData(wb, "PL23_Ctrl_vs_R848_Ctrl", PL23_v_R848)
      writeData(wb, "PL23_2DG_vs_R848_2DG", PL23_2DG_v_R848_2DG)
    
      saveWorkbook(wb, filename, overwrite = TRUE)
    }
    
    return(list(PL23_2DG_v_Ctrl = PL23_2DG_v_Ctrl, 
                R848_2DG_v_Ctrl = R848_2DG_v_Ctrl, 
                PL23_v_R848 = PL23_v_R848, 
                PL23_2DG_v_R848_2DG = PL23_2DG_v_R848_2DG))
  }
  
GSEA_dotplot <- function(gsea_dat, plot_title, NES_direction){
  if(NES_direction == "up") tags_high_color <- "red"
  else if(NES_direction == "down") tags_high_color <- "blue"
  else if(NES_direction == "both") tags_high_color <- "black"
  else {
    print("error: NES_direction must be equal to up, down, or both")
    return(NULL)
  }
  
  if(max(gsea_dat$FDR.q.val) <= 0.05) x_max <- 0.05
  else x_max <- max(gsea_dat$FDR.q.val)
  
  ggplot(gsea_dat,
         aes(x = FDR.q.val, y = reorder(NAME, FDR.q.val, decreasing = TRUE),
             size = SIZE, color = TAGS)) +
    geom_point() +
    scale_size_area(max_size = 5, limits = c(1, max(gsea_dat$SIZE))) +
    scale_colour_gradient(name = "Tags (%)", 
                          low = "white", high = tags_high_color, 
                          limits = c(0, max(gsea_dat$TAGS))) +
    xlab("FDR") +
    scale_x_continuous(limits = c(0.00, x_max)) +
    ylab("") +
    ggtitle(plot_title) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
}

GSEA_dotplot2 <- function(gsea_dat, plot_title, name_type, NES_min, NES_max){
  if(nrow(gsea_dat) > 35){
    gsea_dat <- gsea_dat[order(gsea_dat$FDR.q.val, decreasing = FALSE), ]
    gsea_dat <- gsea_dat[1:35, ]
    plot_title <- stri_join(c(plot_title, "(Top 35)"), collapse = "\n")
  }
  
  if(max(gsea_dat$FDR.q.val) <= 0.01) x_max <- max(gsea_dat$FDR.q.val)
  else if(max(gsea_dat$FDR.q.val) <= 0.05) x_max <- 0.05
  else x_max <- max(gsea_dat$FDR.q.val)
  
  if(missing(name_type)) name_type <- "na"
  else if(name_type == "full") gsea_dat$NAME <- gsea_dat$FULL_NAME
  
  ggplot(gsea_dat,
         aes(x = FDR.q.val, y = reorder(NAME, FDR.q.val, decreasing = TRUE),
             size = SIZE, color = NES)) +
    geom_point() +
    scale_size_area(max_size = 5, limits = c(1, max(gsea_dat$SIZE))) +
    scale_colour_gradient2(name = "NES", 
                          low = "blue", mid = "white", high = "red", 
                          limits = c(NES_min, NES_max)) +
    xlab("FDR") +
    scale_x_continuous(limits = c(0.00, x_max)) +
    ylab("") +
    ggtitle(plot_title) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
}



# Hallmark Canonical Pathways --------------------------------------------------
  # hallmark <- compile_results("handling_GSEA_results/test.xlsx",
  hallmark <- compile_results(NULL,
                              PL23_2DG_v_Ctrl_up = "handling_GSEA_results/GSEA_results/Hallmark_PL23_2DG_vs_PL23/gsea_report_for_PL2-3_2DG.tsv",
                              PL23_2DG_v_Ctrl_down = "handling_GSEA_results/GSEA_results/Hallmark_PL23_2DG_vs_PL23/gsea_report_for_PL2-3.tsv",
                              R848_2DG_v_Ctrl_up = "handling_GSEA_results/GSEA_results/Hallmark_R848_2DG_vs_R848/gsea_report_for_R848_2DG.tsv",
                              R848_2DG_v_Ctrl_down = "handling_GSEA_results/GSEA_results/Hallmark_R848_2DG_vs_R848/gsea_report_for_R848.tsv",
                              PL23_v_R848_up = "handling_GSEA_results/GSEA_results/Hallmark_PL23_vs_R848/gsea_report_for_PL2-3.tsv",
                              PL23_v_R848_down = "handling_GSEA_results/GSEA_results/Hallmark_PL23_vs_R848/gsea_report_for_R848.tsv",
                              PL23_2DG_v_R848_2DG_up = "handling_GSEA_results/GSEA_results/Hallmark_PL23_2DG_vs_R848_2DG/gsea_report_for_PL2-3_2DG.tsv",
                              PL23_2DG_v_R848_2DG_down = "handling_GSEA_results/GSEA_results/Hallmark_PL23_2DG_vs_R848_2DG/gsea_report_for_R848_2DG.tsv")
  head(hallmark$PL23_2DG_v_Ctrl)
  

  # Hallmark pathways dotplots -------------------------------------------------
    # hallmark_up <- GSEA_dotplot(subset(hallmark$PL23_2DG_v_Ctrl, NES > 0 & FDR.q.val < 0.05),
    #                             plot_title = "PL2-3: 2DG vs Control\nUpregulated Hallmark Canonical Pathways",
    #                             NES_direction = "up")
    # hallmark_up
    # ggsave("handling_GSEA_results/GSEA_dotplots/Hallmark_PL23_2DG_vs_PL23_up.png", 
    #        plot = hallmark_up, 
    #        units = "in", width = 6, height = 4, dpi = 300)
    # 
    # hallmark_down <- GSEA_dotplot(subset(hallmark$PL23_2DG_v_Ctrl, NES < 0 & FDR.q.val < 0.05),
    #                               plot_title = "PL2-3: 2DG vs Control\nDownregulated Hallmark Canonical Pathways",
    #                               NES_direction = "down")
    # hallmark_down
    # ggsave("handling_GSEA_results/GSEA_dotplots/Hallmark_PL23_2DG_vs_PL23_down.png", 
    #        plot = hallmark_down, 
    #        units = "in", width = 6, height = 8, dpi = 300)
    
    min(subset(hallmark$PL23_2DG_v_Ctrl, FDR.q.val < 0.05)$NES)
    max(subset(hallmark$PL23_2DG_v_Ctrl, FDR.q.val < 0.05)$NES)
    hallmark_PL23 <- GSEA_dotplot2(subset(hallmark$PL23_2DG_v_Ctrl, FDR.q.val < 0.05),
                                   plot_title = "PL2-3: 2DG vs Control\nHallmark Canonical Pathways",
                                   NES_min = -3, NES_max = 3)
    hallmark_PL23
    ggsave("handling_GSEA_results/GSEA_dotplots/Hallmark_PL23_2DG_vs_PL23.png", 
           plot = hallmark_PL23, 
           units = "in", width = 6, height = 8, dpi = 300)
    
    
    min(subset(hallmark$R848_2DG_v_Ctrl, FDR.q.val < 0.05)$NES)
    max(subset(hallmark$R848_2DG_v_Ctrl, FDR.q.val < 0.05)$NES)
    hallmark_R848 <- GSEA_dotplot2(subset(hallmark$R848_2DG_v_Ctrl, FDR.q.val < 0.05),
                                   plot_title = "R848: 2DG vs Control\nHallmark Canonical Pathways",
                                   NES_min = -3, NES_max = 3)
    hallmark_R848
    ggsave("handling_GSEA_results/GSEA_dotplots/Hallmark_R848_2DG_vs_R848.png", 
           plot = hallmark_R848, 
           units = "in", width = 6, height = 6, dpi = 300)
    
    
    min(subset(hallmark$PL23_v_R848, FDR.q.val < 0.05)$NES)
    max(subset(hallmark$PL23_v_R848, FDR.q.val < 0.05)$NES)
    hallmark_PL23_v_R848 <- GSEA_dotplot2(subset(hallmark$PL23_v_R848, FDR.q.val < 0.05),
                                          plot_title = "PL2-3 vs R848\nHallmark Canonical Pathways",
                                          NES_min = -3.2, NES_max = 3)
    hallmark_PL23_v_R848
    ggsave("handling_GSEA_results/GSEA_dotplots/Hallmark_PL23_vs_R848.png", 
           plot = hallmark_PL23_v_R848, 
           units = "in", width = 6, height = 8, dpi = 300)
    
# All Canonical Pathways
  # canonicalPWs <- compile_results("handling_GSEA_results/CanonicalPathways_Results.xlsx",
  canonicalPWs <- compile_results(NULL,
                                  PL23_2DG_v_Ctrl_up = "handling_GSEA_results/GSEA_results/CanonicalPathways_all_PL23_2DG_vs_PL23/gsea_report_for_PL2-3_2DG.tsv",
                                  PL23_2DG_v_Ctrl_down = "handling_GSEA_results/GSEA_results/CanonicalPathways_all_PL23_2DG_vs_PL23/gsea_report_for_PL2-3.tsv",
                                  R848_2DG_v_Ctrl_up = "handling_GSEA_results/GSEA_results/CanonicalPathways_all_R848_2DG_vs_R848/gsea_report_for_R848_2DG.tsv",
                                  R848_2DG_v_Ctrl_down = "handling_GSEA_results/GSEA_results/CanonicalPathways_all_R848_2DG_vs_R848/gsea_report_for_R848.tsv",
                                  PL23_v_R848_up = "handling_GSEA_results/GSEA_results/CanonicalPathways_all_PL23_vs_R848/gsea_report_for_PL2-3.tsv",
                                  PL23_v_R848_down = "handling_GSEA_results/GSEA_results/CanonicalPathways_all_PL23_vs_R848/gsea_report_for_R848.tsv",
                                  PL23_2DG_v_R848_2DG_up = "handling_GSEA_results/GSEA_results/CanonicalPathways_all_PL23_2DG_vs_R848_2DG/gsea_report_for_PL2-3_2DG.tsv",
                                  PL23_2DG_v_R848_2DG_down = "handling_GSEA_results/GSEA_results/CanonicalPathways_all_PL23_2DG_vs_R848_2DG/gsea_report_for_R848_2DG.tsv")
  
  min(subset(canonicalPWs$PL23_2DG_v_Ctrl, FDR.q.val < 0.05)$NES)
  max(subset(canonicalPWs$PL23_2DG_v_Ctrl, FDR.q.val < 0.05)$NES)
  canonicalPWs_PL23 <- GSEA_dotplot2(subset(canonicalPWs$PL23_2DG_v_Ctrl, FDR.q.val < 0.05),
                                     plot_title = "PL2-3: 2DG vs Control\nCanonical Pathways",
                                     name_type = "full", NES_min = -3, NES_max = 3)
  canonicalPWs_PL23
  ggsave("handling_GSEA_results/GSEA_dotplots/canonicalPWs_PL23_2DG_vs_PL23.png", 
         plot = canonicalPWs_PL23, 
         units = "in", width = 12.5, height = 8, dpi = 300)
  
  reactome_PL23 <- GSEA_dotplot2(subset(canonicalPWs$PL23_2DG_v_Ctrl, FDR.q.val < 0.05 & SOURCE == "REACTOME"),
                                 plot_title = "PL2-3: 2DG vs Control\nReactome Canonical Pathways",
                                 NES_min = -3, NES_max = 3)
  reactome_PL23
  ggsave("handling_GSEA_results/GSEA_dotplots/Reactome_PL23_2DG_vs_PL23.png", 
         plot = reactome_PL23, 
         units = "in", width = 11, height = 8, dpi = 300)
  
  
  
  min(subset(canonicalPWs$R848_2DG_v_Ctrl, FDR.q.val < 0.05)$NES)
  max(subset(canonicalPWs$R848_2DG_v_Ctrl, FDR.q.val < 0.05)$NES)
  canonicalPWs_R848 <- GSEA_dotplot2(subset(canonicalPWs$R848_2DG_v_Ctrl, FDR.q.val < 0.05),
                                 plot_title = "R848: 2DG vs Control\nCanonical Pathways",
                                 name_type = "full", NES_min = -3, NES_max = 3)
  canonicalPWs_R848
  ggsave("handling_GSEA_results/GSEA_dotplots/canonicalPWs_R848_2DG_vs_R848.png", 
         plot = canonicalPWs_R848, 
         units = "in", width = 13, height = 8, dpi = 300)
  
  reactome_R848 <- GSEA_dotplot2(subset(canonicalPWs$R848_2DG_v_Ctrl, FDR.q.val < 0.05 & SOURCE == "REACTOME"),
                                 plot_title = "R848: 2DG vs Control\nReactome Canonical Pathways",
                                 NES_min = -3, NES_max = 3)
  reactome_R848
  ggsave("handling_GSEA_results/GSEA_dotplots/Reactome_R848_2DG_vs_R848.png", 
         plot = reactome_R848, 
         units = "in", width = 12.5, height = 8, dpi = 300)
  
  
  min(subset(canonicalPWs$PL23_v_R848, FDR.q.val < 0.05)$NES)
  max(subset(canonicalPWs$PL23_v_R848, FDR.q.val < 0.05)$NES)
  canonicalPWs_PL23_v_R848 <- GSEA_dotplot2(subset(canonicalPWs$PL23_v_R848, FDR.q.val < 0.05),
                                        plot_title = "PL2-3 vs R848\nCanonical Pathways",
                                        name_type = "full", NES_min = -3, NES_max = 3)
  canonicalPWs_PL23_v_R848
  ggsave("handling_GSEA_results/GSEA_dotplots/canonicalPWs_PL23_vs_R848.png", 
         plot = canonicalPWs_PL23_v_R848, 
         units = "in", width = 13, height = 8, dpi = 300)
  
  reactome_PL23_v_R848 <- GSEA_dotplot2(subset(canonicalPWs$PL23_v_R848, FDR.q.val < 0.05 & SOURCE == "REACTOME"),
                                 plot_title = "PL2-3 vs R848\nReactome Canonical Pathways",
                                 NES_min = -3, NES_max = 3)
  reactome_PL23_v_R848
  ggsave("handling_GSEA_results/GSEA_dotplots/Reactome_PL23_vs_R848.png", 
         plot = reactome_PL23_v_R848, 
         units = "in", width = 12.5, height = 8, dpi = 300)
  