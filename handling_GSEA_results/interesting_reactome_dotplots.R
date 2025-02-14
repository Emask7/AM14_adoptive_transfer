reactome_PL23_interesting <- read.xlsx("Output_DO_NOT_OVERWRITE/CanonicalPathways_Results.xlsx",
                                       sheet = "PL23_2DG_vs_Ctrl_reactome_int")
head(reactome_PL23_interesting)


reactome_PL23_int_plot <- GSEA_dotplot2(reactome_PL23_interesting,
                                        plot_title = "PL2-3: 2DG vs Control\nReactome Pathways",
                                        NES_min = -3, NES_max = 3,
                                        ignore_limit = TRUE)
reactome_PL23_int_plot
ggsave("handling_GSEA_results/GSEA_dotplots/Reactome_PL23_2DG_vs_PL23_interesting.png", 
       plot = reactome_PL23_int_plot, 
       units = "in", width = 9, height = 9, dpi = 300)
