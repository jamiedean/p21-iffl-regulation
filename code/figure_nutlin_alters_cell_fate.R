setwd('/Users/jamiedean/Documents/Papers/Code for Publications/p21-iffl-regulation')

source('code/analysis_functions/functional_analysis.R')
source('code/analysis_functions/preprocessing.R')
source('code/analysis_functions/radiosensitivity.R')
source('code/analysis_functions/rna_seq_differential_expression.R')
source('code/analysis_functions/transcription_factor_binding_site_enrichment_analysis.R')
source('code/analysis_functions/visualize_dynamics.R')

###############################################################################################################

# RNA-seq differential expression

# Differential expression
directory <- file.path('data/rna_seq/hafner2017')
csvfile <- file.path(directory, 'sample_table.csv')
coldata <- read.csv(csvfile, stringsAsFactors = FALSE)

coldata$names <- coldata$Run
coldata$files <- file.path(directory, 'quants', paste0(coldata$names, '_quant'), 'quant.sf')

# Run this line if there is no matched transcriptome
#make_linked_transcriptome()

differential_expression_results <-
  perform_differential_expression_analysis(coldata[coldata$time >= 3 & coldata$time <= 12,],
                                           'ir_nutlin_vs_ir_timecourse')

differential_expression_results <- differential_expression_results[complete.cases(differential_expression_results),]
# Number of upregulated genes: 1473
nrow(differential_expression_results[
  differential_expression_results$padj < 0.05 & differential_expression_results$log2FoldChange > 0,])
# Number of downregulated genes: 1188
nrow(differential_expression_results[
  differential_expression_results$padj < 0.05 & differential_expression_results$log2FoldChange < 0,])

p_volcano_rna_seq <- plot_volcano_plot_rna_seq_timecourse(differential_expression_results)

upregulated_genes <- c(t(subset(differential_expression_results, log2FoldChange > 0.5 & padj < 0.05)$symbol))
downregulated_genes <- c(t(subset(differential_expression_results, log2FoldChange < -0.5 & padj < 0.05)$symbol))

write.csv(upregulated_genes, 'data/signatures/ir_nutlin_vs_ir_upregulated_genes.csv')
write.csv(downregulated_genes, 'data/signatures/ir_nutlin_vs_ir_downregulated_genes.csv')

###############################################################################################################

# Transcription factor enrichment analysis

p_tf_enrichment_upregulated_genes <-
  plot_transcription_factor_enrichment(differential_expression_results, 'upregulated')
p_tf_enrichment_downregulated_genes <-
  plot_transcription_factor_enrichment(differential_expression_results, 'downregulated')

###############################################################################################################

# Gene set variation analysis

dds <- create_differential_expression_dataset(coldata[coldata$time >= 3 & coldata$time <= 12,], 'ir_nutlin_vs_ir_timecourse')
dds_norm <- vst(dds)

p_gsva_persister_cells <- plot_geneset_variation_analysis_heatmap(dds_norm, 'custom', 'persister_cells', FALSE, TRUE)

###############################################################################################################

# Signature association with radiation response

p_cell_line_radiotherapy <- signature_outcome_association_cell_lines_ccle()
# Upregulated signature: p = 6.444401e-06
# Downregulated signature: p = 1.455263e-10

metabric_gsva <- compute_geneset_variation_analysis_scores_metabric('ir_nutlin_vs_ir_signature')

p_survival_no_adjuvant_therapy <- 
  signature_outcome_association_clinical_metabric(metabric_gsva, treatment = 'no_adjuvant_therapy')
p_survival_chemo_radiotherapy <- 
  signature_outcome_association_clinical_metabric(metabric_gsva, treatment = 'chemo_radiotherapy')

p_no_adjuvant_therapy_rfs_up <- p_survival_no_adjuvant_therapy[[1]]
p_no_adjuvant_therapy_os_up <- p_survival_no_adjuvant_therapy[[2]]
p_no_adjuvant_therapy_rfs_down <- p_survival_no_adjuvant_therapy[[3]]
p_no_adjuvant_therapy_os_down <- p_survival_no_adjuvant_therapy[[4]]

p_chemo_radiotherapy_rfs_up <- p_survival_chemo_radiotherapy[[1]]
p_chemo_radiotherapy_os_up <- p_survival_chemo_radiotherapy[[2]]
p_chemo_radiotherapy_rfs_down <- p_survival_chemo_radiotherapy[[3]]
p_chemo_radiotherapy_os_down <- p_survival_chemo_radiotherapy[[4]]

###############################################################################################################

# Figure

p_rna_seq <- 
  plot_grid(
    p_volcano_rna_seq,
    p_tf_enrichment_upregulated_genes,
    p_tf_enrichment_downregulated_genes,
    labels = 'AUTO',
    ncol = 3)

p_gsva <- 
  plot_grid(
    grid.grabExpr(draw(p_gsva_persister_cells,
                       heatmap_legend_side = 'left',
                       annotation_legend_side = 'left',
                       padding = unit(c(2, 2, 2, 20), 'mm'))),
    labels = c('D'),
    ncol = 1,
    rel_heights = c(1))

p_radiosensitivity <-
  plot_grid(
    p_cell_line_radiotherapy,
    plot_grid(
      p_chemo_radiotherapy_rfs_down$plot + 
        scale_color_manual(labels = c('< Median', '>= Median'), values = c(plot_colors[1], plot_colors[2])) +
        theme(legend.position = c(0.8, 0.9)),
      p_chemo_radiotherapy_os_down$plot +
        scale_color_manual(labels = c('< Median', '>= Median'), values = c(plot_colors[1], plot_colors[2])) +
        theme(legend.position = c(0.8, 0.9)),
      labels = c('F', 'G')),
    labels = c('E', ''),
    ncol = 1,
    rel_heights = c(1, 1))

figure <- 
  plot_grid(
    p_rna_seq,
    p_gsva,
    p_radiosensitivity,
    ncol = 1,
    rel_heights = c(1, 1, 3))

save_plot('figures/figure_nutlin_alters_cell_fate.pdf', figure, ncol = 2.5, nrow = 5.5)

supplemental_figure <- 
  plot_grid(
    grid.grabExpr(draw(p_gsva_hallmark,
                       heatmap_legend_side = 'left',
                       annotation_legend_side = 'left',
                       padding = unit(c(2, 2, 2, 20), 'mm'))),
    plot_grid(
      p_no_adjuvant_therapy_rfs_down$plot + 
        scale_color_manual(labels = c('< Median', '>= Median'), values = c(plot_colors[1], plot_colors[2])) +
        theme(legend.position = c(0.8, 0.9)),
      p_no_adjuvant_therapy_os_down$plot + 
        scale_color_manual(labels = c('< Median', '>= Median'), values = c(plot_colors[1], plot_colors[2])) +
        theme(legend.position = c(0.8, 0.9)),
      p_no_adjuvant_therapy_rfs_up$plot + 
        scale_color_manual(labels = c('< Median', '>= Median'), values = c(plot_colors[1], plot_colors[2])) +
        theme(legend.position = c(0.8, 0.9)),
      p_no_adjuvant_therapy_os_up$plot + 
        scale_color_manual(labels = c('< Median', '>= Median'), values = c(plot_colors[1], plot_colors[2])) +
        theme(legend.position = c(0.8, 0.9)),
      p_chemo_radiotherapy_rfs_up$plot + 
        scale_color_manual(labels = c('< Median', '>= Median'), values = c(plot_colors[1], plot_colors[2])) +
        theme(legend.position = c(0.8, 0.9)),
      p_chemo_radiotherapy_os_up$plot + 
        scale_color_manual(labels = c('< Median', '>= Median'), values = c(plot_colors[1], plot_colors[2])) +
        theme(legend.position = c(0.8, 0.9)),
      labels = c('B', 'C', 'D', 'E', 'F', 'G'),
      ncol = 2),
    labels = c('A', ''),
    ncol = 1,
    rel_heights = c(3, 3))

save_plot('figures/supplemental_figure_nutlin_alters_cell_fate.pdf', supplemental_figure, ncol = 3, nrow = 8)
