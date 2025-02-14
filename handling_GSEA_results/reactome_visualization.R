# SLC_MEDIATED_TRANSMEMBRANE_TRANSPORT
# COMPLEMENT_CASCADE
# CLASSICAL_ANTIBODY_MEDIATED_COMPLEMENT_ACTIVATION
# INITIAL_TRIGGERING_OF_COMPLEMENT
# COMPLEX_I_BIOGENESIS
# RESPIRATORY_ELECTRON_TRANSPORT
# ESTROGEN_DEPENDENT_GENE_EXPRESSION
# MITOCHONDRIAL_BIOGENESIS
# AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT
# MITOCHONDRIAL_TRANSLATION
# COMPLEX_IV_ASSEMBLY
# METABOLISM_OF_LIPIDS
# IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL
# AMINO_ACID_TRANSPORT_ACROSS_THE_PLASMA_MEMBRANE

save_pw_plot <- function(pathway, file_name, dat, w, h, return_plot){
  dat <- subset(dat, padj <= 0.05)
  gene_list <- dat$log2FoldChange
  names(gene_list) <- dat$entrezgene
  gene_list <- subset(gene_list, !is.na(names(gene_list)))
  dat <- sort(gene_list, decreasing = TRUE)
  pw_plot <- viewPathway(pathway, organism = "mouse", foldChange = dat)
  
  if(missing(w)) w <- 3000
  if(missing(h)) h <- 3000
  
  ggsave(stri_join(c("handling_GSEA_results/Pathway_Diagrams/", file_name),
                   collapse = ""), 
         units = "px", width = w, height = h, dpi = 300)
  if(missing(return_plot)) return_plot <- FALSE
  else if(return_plot != TRUE & return_plot != FALSE) return_plot <- FALSE
  
  if(return_plot) return(pw_plot)
}



save_pw_plot("Complement cascade", 
             "REACTOME_COMPLEMENT_CASCADE_PL23_2DGvCtrl.png", 
             res_PL23_full)
save_pw_plot("Classical antibody-mediated complement activation", 
             "REACTOME_CLASSICAL_ANTIBODY_MEDIATED_COMPLEMENT_ACTIVATION_PL23_2DGvCtrl.png", 
             res_PL23_full)
save_pw_plot("Initial triggering of complement", 
             "REACTOME_INITIAL_TRIGGERING_OF_COMPLEMENT_PL23_2DGvCtrl.png", 
             res_PL23_full)
save_pw_plot("Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell", 
             "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL_PL23_2DGvCtrl.png", 
             res_PL23_full)

save_pw_plot("Complex I biogenesis", 
             "REACTOME_COMPLEX_I_BIOGENESIS_PL23_2DGvCtrl.png", 
             res_PL23_full)
save_pw_plot("Respiratory electron transport", 
             "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT_PL23_2DGvCtrl.png", 
             res_PL23_full)
save_pw_plot("Mitochondrial biogenesis", 
             "REACTOME_MITOCHONDRIAL_BIOGENESIS_PL23_2DGvCtrl.png", 
             res_PL23_full)
save_pw_plot("Aerobic respiration and respiratory electron transport", 
             "REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT_PL23_2DGvCtrl.png", 
             res_PL23_full)
save_pw_plot("Mitochondrial translation", 
             "REACTOME_MITOCHONDRIAL_TRANSLATION_PL23_2DGvCtrl.png", 
             res_PL23_full)
save_pw_plot("Complex IV assembly", 
             "REACTOME_COMPLEX_IV_ASSEMBLY_PL23_2DGvCtrl.png", 
             res_PL23_full)

save_pw_plot("SLC-mediated transmembrane transport", 
             "REACTOME_SLC_MEDIATED_TRANSMEMBRANE_TRANSPORT_PL23_2DGvCtrl.png", 
             res_PL23_full, w = 4000, h = 4000, return_plot = FALSE)
save_pw_plot("Amino acid transport across the plasma membrane", 
             "REACTOME_AMINO_ACID_TRANSPORT_ACROSS_THE_PLASMA_MEMBRANE_PL23_2DGvCtrl.png", 
             res_PL23_full)

# save_pw_plot("Estrogen-dependent gene expression", 
#              "REACTOME_ESTROGEN_DEPENDENT_GENE_EXPRESSION_PL23_2DGvCtrl.png", 
#              res_PL23_full)
save_pw_plot("Metabolism of lipids", 
             "REACTOME_METABOLISM_OF_LIPIDS_PL23_2DGvCtrl.png", 
             res_PL23_full)







save_pw_plot("Complement cascade", 
             "REACTOME_COMPLEMENT_CASCADE_R848_2DGvCtrl.png", 
             res_R848_full)
save_pw_plot("Classical antibody-mediated complement activation", 
             "REACTOME_CLASSICAL_ANTIBODY_MEDIATED_COMPLEMENT_ACTIVATION_R848_2DGvCtrl.png", 
             res_R848_full)
save_pw_plot("Initial triggering of complement", 
             "REACTOME_INITIAL_TRIGGERING_OF_COMPLEMENT_R848_2DGvCtrl.png", 
             res_R848_full)
save_pw_plot("Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell", 
             "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL_R848_2DGvCtrl.png", 
             res_R848_full)

save_pw_plot("SLC-mediated transmembrane transport", 
             "REACTOME_SLC_MEDIATED_TRANSMEMBRANE_TRANSPORT_R848_2DGvCtrl.png", 
             res_R848_full, w = 4000, h = 4000, return_plot = FALSE)





save_pw_plot("Complement cascade", 
             "REACTOME_COMPLEMENT_CASCADE_PL23vR848.png", 
             res_PL23vR848_full)
save_pw_plot("Classical antibody-mediated complement activation", 
             "REACTOME_CLASSICAL_ANTIBODY_MEDIATED_COMPLEMENT_ACTIVATION_PL23vR848.png", 
             res_PL23vR848_full)
save_pw_plot("Initial triggering of complement", 
             "REACTOME_INITIAL_TRIGGERING_OF_COMPLEMENT_PL23vR848.png", 
             res_PL23vR848_full)
save_pw_plot("Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell", 
             "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL_PL23vR848.png", 
             res_PL23vR848_full)

save_pw_plot("SLC-mediated transmembrane transport", 
             "REACTOME_SLC_MEDIATED_TRANSMEMBRANE_TRANSPORT_PL23vR848.png", 
             res_PL23vR848_full, w = 4000, h = 4000, return_plot = FALSE)
