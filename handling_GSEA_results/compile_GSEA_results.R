# Function that takes TSV results from all 4 comparisons and compiles them into 1 Excel file
  compile_results <- function(filename, 
                              PL23_2DG_v_Ctrl_up, R848_2DG_v_Ctrl_up, PL23_v_R848_up, PL23_2DG_v_R848_2DG_up, 
                              PL23_2DG_v_Ctrl_down, R848_2DG_v_Ctrl_down, PL23_v_R848_down, PL23_2DG_v_R848_2DG_down){
    PL23_2DG_v_Ctrl_up <- read.delim(PL23_2DG_v_Ctrl_up)
    R848_2DG_v_Ctrl_up <- read.delim(R848_2DG_v_Ctrl_up) 
    PL23_v_R848_up <- read.delim(PL23_v_R848_up)
    PL23_2DG_v_R848_2DG_up <- read.delim(PL23_2DG_v_R848_2DG_up)
    PL23_2DG_v_Ctrl_down <- read.delim(PL23_2DG_v_Ctrl_down)
    R848_2DG_v_Ctrl_down <- read.delim(R848_2DG_v_Ctrl_down)
    PL23_v_R848_down <- read.delim(PL23_v_R848_down)
    PL23_2DG_v_R848_2DG_down <- read.delim(PL23_2DG_v_R848_2DG_down)
    
    PL23_2DG_v_Ctrl_all <- full_join(PL23_2DG_v_Ctrl_up, PL23_2DG_v_Ctrl_down)
    PL23_2DG_v_Ctrl_full <- PL23_2DG_v_Ctrl_full[, c(1, 4:11)]
    
    R848_2DG_v_Ctrl_all <- full_join(R848_2DG_v_Ctrl_up, R848_2DG_v_Ctrl_down)
    R848_2DG_v_Ctrl_all <- R848_2DG_v_Ctrl_all[, c(1, 4:11)]
    
    PL23_v_R848_all <- full_join(PL23_v_R848_up, PL23_v_R848_down)
    PL23_v_R848_all <- PL23_v_R848_all[, c(1, 4:11)]
    
    PL23_2DG_v_R848_2DG_all <- full_join(PL23_2DG_v_R848_2DG_up, PL23_2DG_v_R848_2DG_down)
    PL23_2DG_v_R848_2DG_all <- PL23_2DG_v_R848_2DG_all[, c(1, 4:11)]
    
    wb <- createWorkbook(filename)
    
    addWorksheet(wb, "PL23_2DG_vs_Ctrl")
    addWorksheet(wb, "R848_2DG_vs_Ctrl")
    addWorksheet(wb, "PL23_Ctrl_vs_R848_Ctrl")
    addWorksheet(wb, "PL23_2DG_vs_R848_2DG")
    
    writeData(wb, "PL23_2DG_vs_Ctrl", PL23_2DG_v_Ctrl_all)
    writeData(wb, "R848_2DG_vs_Ctrl", R848_2DG_v_Ctrl_all)
    writeData(wb, "PL23_Ctrl_vs_R848_Ctrl", PL23_v_R848_all)
    writeData(wb, "PL23_2DG_vs_R848_2DG", PL23_2DG_v_R848_2DG_all)
    
    saveWorkbook(wb, filename, overwrite = TRUE)
  }
  
  # Hallmark Canonical Pathways
  compile_results("handling_GSEA_results/Hallmark_Results.xlsx"
                  PL23_2DG_v_Ctrl_up = read.delim(),
                  R848_2DG_v_Ctrl_up = read.delim(), 
                  PL23_v_R848_up = read.delim(),
                  PL23_2DG_v_R848_2DG_up = read.delim(),
                  PL23_2DG_v_Ctrl_down = read.delim(),
                  R848_2DG_v_Ctrl_down = read.delim(),
                  PL23_v_R848_down = read.delim(),
                  PL23_2DG_v_R848_2DG_down = read.delim(),
  )
  
  
  
  
  
  PL23_2DG_v_Ctrl_up <- read.delim("handling_GSEA_results/GSEA_results/Hallmark_PL23_2DG_vs_PL23/gsea_report_for_PL2-3_2DG_1737491393564.tsv")
  PL23_2DG_v_Ctrl_down <- read.delim("handling_GSEA_results/GSEA_results/Hallmark_PL23_2DG_vs_PL23/gsea_report_for_PL2-3_1737491393564.tsv")
  
  nrow(PL23_2DG_v_Ctrl_up) #8
  nrow(PL23_2DG_v_Ctrl_down) #42
  
  PL23_2DG_v_Ctrl_full <- full_join(PL23_2DG_v_Ctrl_up, PL23_2DG_v_Ctrl_down)
  nrow(PL23_2DG_v_Ctrl_full)
  head(PL23_2DG_v_Ctrl_full)
  PL23_2DG_v_Ctrl_full <- PL23_2DG_v_Ctrl_full[, c(1, 4:11)]
  
  