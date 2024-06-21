setwd('/Users/jamiedean/Documents/Papers/Code for Publications/p21-iffl-regulation')

source('code/analysis_functions/preprocessing.R')
source('code/analysis_functions/visualize_dynamics.R')

################################################################################################################

experiment <- 'IR_Nutlin'

p53 <- load_data(experiment, 'traces_CFP', FALSE)
p21rna <- load_data(experiment, 'traces_YFP', FALSE)
p21 <- load_data(experiment, 'traces_Texas', FALSE)

p53_heatmap <- plot_single_cell_dynamics_heatmaps(p53, 'p53')
png('figures/ir_nutlin_p53_heatmap.png', width = 9, height = 10, units = 'in', res = 1200)
p53_heatmap
dev.off()
p_p53_heatmap <- ggdraw() + 
  draw_image('figures/ir_nutlin_p53_heatmap.png', x = 1, width = 1, height = 1, hjust = 1)

p21ms2_heatmap <- plot_single_cell_dynamics_heatmaps(p21rna, 'p21-MS2')
png('figures/ir_nutlin_p21ms2_heatmap.png', width = 9, height = 10, units = 'in', res = 1200)
p21ms2_heatmap
dev.off()
p_p21ms2_heatmap <- ggdraw() + 
  draw_image('figures/ir_nutlin_p21ms2_heatmap.png', x = 1, width = 1, height = 1, hjust = 1)

p21_heatmap <- plot_single_cell_dynamics_heatmaps(p21, 'p21')
png('figures/ir_nutlin_p21_heatmap.png', width = 9, height = 10, units = 'in', res = 1200)
p21_heatmap
dev.off()
p_p21_heatmap <- ggdraw() + 
  draw_image('figures/ir_nutlin_p21_heatmap.png', x = 1, width = 1, height = 1, hjust = 1)

################################################################################################################

figure <- 
  plot_grid(
    p_p53_heatmap,
    p_p21ms2_heatmap,
    p_p21_heatmap,
    ncol = 3,
    labels = 'AUTO'
  )

save_plot('figures/supplemental_figure_model_validation.pdf', figure,
          ncol = 3, nrow = 2)
