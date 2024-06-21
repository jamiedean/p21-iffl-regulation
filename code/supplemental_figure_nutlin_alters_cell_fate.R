setwd('/Users/jamiedean/Documents/Papers/Code for Publications/p21-iffl-regulation')

source('code/analysis_functions/functional_analysis.R')
source('code/analysis_functions/radiosensitivity.R')
source('code/analysis_functions/rna_seq_differential_expression.R')

###############################################################################################################

# Gene set variation analysis

directory <- file.path('data/rna_seq/hafner2017')
csvfile <- file.path(directory, 'sample_table.csv')
coldata <- read.csv(csvfile, stringsAsFactors = FALSE)

coldata$names <- coldata$Run
coldata$files <- file.path(directory, 'quants', paste0(coldata$names, '_quant'), 'quant.sf')

# Run this line if there is no matched transcriptome
#make_linked_transcriptome()

dds <- create_differential_expression_dataset(coldata[coldata$time >= 3 & coldata$time <= 12,], 'ir_nutlin_vs_ir_timecourse')
dds_norm <- vst(dds)

p_gsva_hallmark <- plot_geneset_variation_analysis_heatmap(dds_norm, 'custom', 'hallmark', FALSE, TRUE)

###############################################################################################################

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

###############################################################################################################

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
