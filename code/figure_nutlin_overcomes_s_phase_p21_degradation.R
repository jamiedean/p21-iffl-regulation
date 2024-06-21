setwd('/Users/jamiedean/Documents/Papers/Code for Publications/p21-iffl-regulation')

source('code/analysis_functions/functional_analysis.R')
source('code/analysis_functions/infer_production_and_degradation_rates.R')
source('code/analysis_functions/rna_seq_differential_expression.R')

################################################################################################################

experiment <- 'IR_p53siRNA_p21siRNA'

p21 <- load_data(experiment, 'traces_Texas', FALSE)

p_p21_degradation_single_cell <- 
  infer_protein_degradation_rates_single_cell(p21, cell_number = 12)

p21_degradation_rate_distribution <- 
  infer_protein_degradation_rates(p21[p21$condition == 'p53 siRNA p21 siRNA' | p21$condition == 'p21 siRNA',])
p_p21_degradation_rate_distribution <- p21_degradation_rate_distribution[[1]]
p21_degradation_rate_mixture_model_fit <- p21_degradation_rate_distribution[[2]]

production_rate_plots_ir_nutlin <- 
  infer_production_rates('IR_Nutlin', p21_degradation_rate_mixture_model_fit, time_from_end_of_experiment = 2.5)

################################################################################################################

# RNA-seq and geneset variation analysis

# Differential expression
directory <- file.path('data/rna_seq/hafner2017')
csvfile <- file.path(directory, 'sample_table.csv')
coldata <- read.csv(csvfile, stringsAsFactors = FALSE)

coldata$names <- coldata$Run
coldata$files <- file.path(directory, 'quants', paste0(coldata$names, '_quant'), 'quant.sf')

# Run this line if there is no matched transcriptome
#make_linked_transcriptome()

dds <- 
  create_differential_expression_dataset(coldata[coldata$time >= 3 & coldata$time <= 12,], 'ir_nutlin_vs_ir_timecourse')
dds_norm <- vst(dds)

p_gsva_cell_cycle <- plot_geneset_variation_analysis_heatmap(dds_norm, 'custom', 'cell_cycle_xue2020', TRUE, TRUE)

################################################################################################################

# Figure

figure <- 
  plot_grid(
    plot_grid(
      production_rate_plots_ir_nutlin[[1]],
      production_rate_plots_ir_nutlin[[2]],
      plot_grid(p_p21_degradation_single_cell,
                p_p21_degradation_rate_distribution,
                ncol = 2,
                labels = c('C', 'D')),
      labels = c('A', 'B', ''),
      ncol = 1,
      rel_heights = c(1.75, 1, 1)),
    plot_grid(
      production_rate_plots_ir_nutlin[[3]],
      production_rate_plots_ir_nutlin[[4]],
      grid.grabExpr(draw(p_gsva_cell_cycle,
                         heatmap_legend_side = 'left',
                         annotation_legend_side = 'left')),
      labels = c('E', 'F', 'G'),
      ncol = 1,
      rel_heights = c(1, 1, 1)),
      ncol = 2
  )

save_plot('figures/figure_nutlin_overcomes_s_phase_p21_degradation.pdf', figure, ncol = 3.75, nrow = 4)
