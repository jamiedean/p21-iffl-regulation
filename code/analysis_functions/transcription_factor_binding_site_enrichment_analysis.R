# Predict which transcription factors are responsible for the differential expression of genes

library(TFEA.ChIP)
# Internal database contains 1060 datasets from the ENCODE project
# Additional databases are available from https://github.com/LauraPS1/TFEA.ChIP_downloads
load('data/chip_seq/ReMap2022+EnsTSS+CellTypeEnh.Rdata')
set_user_data(MetaData, ChIPDB)

source('code/analysis_functions/custom_plot_theme.R')

################################################################################################################

plot_transcription_factor_enrichment <- function(differential_expression_results, selected_name) {
  
  if (selected_name == 'upregulated') {
    differentially_expressed_genes_selected <- 
      c(t(subset(differential_expression_results, log2FoldChange > 0 & padj < 0.05)$symbol))  
  } else if (selected_name == 'downregulated') {
    differentially_expressed_genes_selected <- 
      c(t(subset(differential_expression_results, log2FoldChange < 0 & padj < 0.05)$symbol)) 
  }
  
  # Extract vector with names of non-responsive genes
  control_genes <- c(t(differential_expression_results[differential_expression_results$padj > 0.5, 'symbol']))
  
  contingency_matrix_list <- contingency_matrix(differentially_expressed_genes_selected, control_genes)
  
  # Generates list of p-values and ORs
  pval_mat <- getCMstats(contingency_matrix_list)
  
  if (selected_name == 'upregulated') {
    highlight <- c('TP53', 'TP63')
    names(highlight) <- c('TP53', 'TP63')
    col <- c(plot_colors[1], plot_colors[3])
  } else if (selected_name == 'downregulated') {
    highlight <- c('E2F4')
    names(highlight) <- c('E2F4')
    col <- c(plot_colors[2])
    #highlight <- c('E2F4', 'LIN9')
    #names(highlight) <- c('E2F4', 'LIN9')
    #col <- c(plot_colors[2], plot_colors[5])
  }
  
  # Plot p-values against ORs highlighting indicated TFs
  plot_CM(pval_mat, specialTF = highlight, TF_colors = col)
  
  p_tf_enrichment <- ggplot() +
    xlab(expression(log[2](OR))) + ylab(expression(-log[10](p[adj]))) +
    theme(legend.position = c(0.2, 0.7))
  
  if (selected_name == 'upregulated') {
    
    p_tf_enrichment <- p_tf_enrichment +
      geom_point(data = pval_mat[pval_mat$TF != 'TP53' & pval_mat$TF != 'TP63',],
                 aes(log2.OR, log10.adj.pVal, color = 'gray')) +
      geom_point(data = pval_mat[pval_mat$TF == 'TP53',],
                 aes(log2.OR, log10.adj.pVal, color = plot_colors[1])) +
      geom_point(data = pval_mat[pval_mat$TF == 'TP63',],
                 aes(log2.OR, log10.adj.pVal, color = plot_colors[3])) +
      scale_color_manual(labels = c('TP53', 'TP63', 'Other'), values = c(plot_colors[1], plot_colors[3], 'gray'))
    
  }  else if (selected_name == 'downregulated') {
    
    p_tf_enrichment <- p_tf_enrichment +
      geom_point(data = pval_mat,#[pval_mat$TF != 'E2F4' & pval_mat$TF != 'LIN9',],
                 aes(log2.OR, log10.adj.pVal, color = 'gray')) +
      geom_point(data = pval_mat[pval_mat$TF == 'E2F4',],
                 aes(log2.OR, log10.adj.pVal, color = plot_colors[2])) +
      #geom_point(data = pval_mat[pval_mat$TF == 'LIN9',],
      #           aes(log2.OR, log10.adj.pVal, color = plot_colors[5])) +
      scale_color_manual(labels = c('E2F4', 'Other'), values = c(plot_colors[2], 'gray'))
      #scale_color_manual(labels = c('E2F4', 'LIN9', 'Other'), values = c(plot_colors[2], plot_colors[5], 'gray'))
  }
  
  return(p_tf_enrichment)
}
