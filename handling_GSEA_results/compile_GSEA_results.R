# Function that takes TSV results from all 4 comparisons and compiles them into 1 Excel file
  combine_results <- function(up, down){
    up <- read.delim(up)
    down <- read.delim(down)
    combined <- full_join(up, down)
    combined <- combined[, c(1, 4:11)]
    combined <- separate_wider_delim(combined, cols = NAME, 
                                     delim = stringr::regex("_"),
                                     names = c("SOURCE", "NAME"), too_many = "merge")

    combined <- separate_wider_delim(combined, 
                                     cols = LEADING.EDGE, 
                                     delim = stringr::regex("%|="),
                                     names = c(NA, "TAGS", NA, "LIST", NA, "SIGNAL", NA))

    return(combined[order(combined$FDR.q.val), ])
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
  
  # Hallmark Canonical Pathways
  compile_results("handling_GSEA_results/Hallmark_Results.xlsx",
                  PL23_2DG_v_Ctrl_up = "handling_GSEA_results/GSEA_results/Hallmark_PL23_2DG_vs_PL23/gsea_report_for_PL2-3_2DG.tsv",
                  PL23_2DG_v_Ctrl_down = "handling_GSEA_results/GSEA_results/Hallmark_PL23_2DG_vs_PL23/gsea_report_for_PL2-3.tsv",
                  R848_2DG_v_Ctrl_up = "handling_GSEA_results/GSEA_results/Hallmark_R848_2DG_vs_R848/gsea_report_for_R848_2DG.tsv",
                  R848_2DG_v_Ctrl_down = "handling_GSEA_results/GSEA_results/Hallmark_R848_2DG_vs_R848/gsea_report_for_R848.tsv",
                  PL23_v_R848_up = "handling_GSEA_results/GSEA_results/Hallmark_PL23_vs_R848/gsea_report_for_PL2-3.tsv",
                  PL23_v_R848_down = "handling_GSEA_results/GSEA_results/Hallmark_PL23_vs_R848/gsea_report_for_R848.tsv",
                  PL23_2DG_v_R848_2DG_up = "handling_GSEA_results/GSEA_results/Hallmark_PL23_2DG_vs_R848_2DG/gsea_report_for_PL2-3_2DG.tsv",
                  PL23_2DG_v_R848_2DG_down = "handling_GSEA_results/GSEA_results/Hallmark_PL23_2DG_vs_R848_2DG/gsea_report_for_R848_2DG.tsv"
  )
  
  
# All Canonical Pathways
  compile_results("handling_GSEA_results/CanonicalPathways_Results.xlsx",
                  PL23_2DG_v_Ctrl_up = "handling_GSEA_results/GSEA_results/CanonicalPathways_all_PL23_2DG_vs_PL23/gsea_report_for_PL2-3_2DG.tsv",
                  PL23_2DG_v_Ctrl_down = "handling_GSEA_results/GSEA_results/CanonicalPathways_all_PL23_2DG_vs_PL23/gsea_report_for_PL2-3.tsv",
                  R848_2DG_v_Ctrl_up = "handling_GSEA_results/GSEA_results/CanonicalPathways_all_R848_2DG_vs_R848/gsea_report_for_R848_2DG.tsv",
                  R848_2DG_v_Ctrl_down = "handling_GSEA_results/GSEA_results/CanonicalPathways_all_R848_2DG_vs_R848/gsea_report_for_R848.tsv",
                  PL23_v_R848_up = "handling_GSEA_results/GSEA_results/CanonicalPathways_all_PL23_vs_R848/gsea_report_for_PL2-3.tsv",
                  PL23_v_R848_down = "handling_GSEA_results/GSEA_results/CanonicalPathways_all_PL23_vs_R848/gsea_report_for_R848.tsv",
                  PL23_2DG_v_R848_2DG_up = "handling_GSEA_results/GSEA_results/CanonicalPathways_all_PL23_2DG_vs_R848_2DG/gsea_report_for_PL2-3_2DG.tsv",
                  PL23_2DG_v_R848_2DG_down = "handling_GSEA_results/GSEA_results/CanonicalPathways_all_PL23_2DG_vs_R848_2DG/gsea_report_for_R848_2DG.tsv"
  )
  
  
