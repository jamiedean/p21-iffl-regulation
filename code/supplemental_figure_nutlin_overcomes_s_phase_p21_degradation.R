setwd('/Users/jamiedean/Documents/Papers/Code for Publications/p21-iffl-regulation')

source('code/analysis_functions/infer_production_and_degradation_rates.R')
source('code/analysis_functions/functional_analysis.R')
source('code/analysis_functions/rna_seq_differential_expression.R')
source('code/analysis_functions/visualize_dynamics.R')

################################################################################################################

experiment <- 'IR_p53siRNA_p21siRNA'

p21 <- load_data(experiment, 'traces_Texas', FALSE)

p21_heatmap <- plot_single_cell_dynamics_heatmaps(
  p21[p21$condition == 'p21 siRNA' | p21$condition == 'p53 siRNA p21 siRNA',], 'p21', column_split = TRUE)
png('figures/p21_heatmap_IR_p53siRNA_p21siRNA.png', width = 9, height = 10, units = 'in', res = 1200)
p21_heatmap
dev.off()
p_p21_heatmap <- ggdraw() + 
  draw_image('figures/p21_heatmap_IR_p53siRNA_p21siRNA.png', x = 1, width = 1, height = 1, hjust = 1)
p_p21_heatmap

p21_degradation_rate_distribution <- 
  infer_protein_degradation_rates(p21[p21$condition == 'p53 siRNA p21 siRNA' | p21$condition == 'p21 siRNA',])
p21_degradation_rate_mixture_model_fit <- p21_degradation_rate_distribution[[2]]

production_rate_plots_ir_nutlin <- 
  infer_production_rates('IR_Nutlin', p21_degradation_rate_mixture_model_fit, time_from_end_of_experiment = 2.5)

###############################################################################################################

# Nuclear size

experiment <- 'IR_Nutlin'

area <- load_data(experiment, 'traces_area', FALSE)

p_area_dynamics <- plot_summary_dynamics_by_condition(area, experiment, return_plot = TRUE)
p_area_dynamics <- p_area_dynamics + ylab('Area (a.u.)')
p_area_change <- plot_change_in_traces(area)

###############################################################################################################

# RNA-seq and geneset variation analysis

directory <- file.path('data/rna_seq/hafner2017')
csvfile <- file.path(directory, 'sample_table.csv')
coldata <- read.csv(csvfile, stringsAsFactors = FALSE)

coldata$names <- coldata$Run
coldata$files <- file.path(directory, 'quants', paste0(coldata$names, '_quant'), 'quant.sf')

# Run this line if there is no matched transcriptome
#make_linked_transcriptome()

dds <- create_differential_expression_dataset(coldata[coldata$time >= 3 & coldata$time <= 12,], 
                                              'ir_nutlin_vs_ir_timecourse')
dds_norm <- vst(dds)

p_gsva_nuclear_size <- plot_geneset_variation_analysis_heatmap(dds_norm, 'custom', 'nuclear_size', FALSE, TRUE)

###############################################################################################################

p_microscopy <- 
  plot_grid(
    p_area_dynamics, p_area_change,
    labels = c('C', 'D'),
    ncol = 2)

supplemental_figure <-
  plot_grid(
    p_p21_heatmap,
    production_rate_plots_ir_nutlin[[5]],
    p_microscopy,
    grid.grabExpr(draw(p_gsva_nuclear_size,
                       heatmap_legend_side = 'left',
                       annotation_legend_side = 'left')),
    labels = c('A', 'B', '', 'E'),
    ncol = 1, rel_heights = c(2, 1.75, 1, 1)
  )

save_plot('figures/supplemental_figure_nutlin_overcomes_s_phase_p21_degradation.pdf',
          supplemental_figure, ncol = 2.5, nrow = 6)
