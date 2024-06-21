setwd('/Users/jamiedean/Documents/Papers/Code for Publications/p21-iffl-regulation')

source('code/analysis_functions/prediction_of_transcription_by_transcription_factor_metrics.R')
source('code/analysis_functions/preprocessing.R')
source('code/analysis_functions/visualize_dynamics.R')

################################################################################################################

experiment <- 'IR_15min_resolution'

p53 <- load_data(experiment, 'traces_CFP', FALSE)
p21rna <- load_data(experiment, 'traces_MS2', FALSE)
p21 <- load_data(experiment, 'traces_Texas', FALSE)

p53_heatmap <- plot_single_cell_dynamics_heatmaps(p53, 'p53')
png('figures/p53_heatmap.png', width = 9, height = 10, units = 'in', res = 1200)
p53_heatmap
dev.off()
p_p53_heatmap <- ggdraw() + 
  draw_image('figures/p53_heatmap.png', x = 1, width = 1, height = 1, hjust = 1)

p21rna_heatmap <- plot_single_cell_dynamics_heatmaps(p21rna, 'p21-MS2')
png('figures/p21rna_heatmap.png', width = 9, height = 10, units = 'in', res = 1200)
p21rna_heatmap
dev.off()
p_p21rna_heatmap <- ggdraw() + 
  draw_image('figures/p21rna_heatmap.png', x = 1, width = 1, height = 1, hjust = 1)

p21_heatmap <- plot_single_cell_dynamics_heatmaps(p21, 'p21')
png('figures/p21_heatmap.png', width = 9, height = 10, units = 'in', res = 1200)
p21_heatmap
dev.off()
p_p21_heatmap <- ggdraw() + 
  draw_image('figures/p21_heatmap.png', x = 1, width = 1, height = 1, hjust = 1)

p_heatmaps <- 
  plot_grid(
    p_p53_heatmap,
    p_p21rna_heatmap,
    p_p21_heatmap,
    ncol = 3,
    labels = 'AUTO')

################################################################################################################

cells <- intersect(unique(p53$variable), unique(p21rna$variable))

example_cell_1_number <- 200
example_cell_2_number <- 36

p53_example_cell_1 <- p53[p53$variable == cells[example_cell_1_number], c('time', 'value')]
p21rna_example_cell_1 <- p21rna[p21rna$variable == cells[example_cell_1_number], c('time', 'value')]
p53_example_cell_2 <- p53[p53$variable == cells[example_cell_2_number], c('time', 'value')]
p21rna_example_cell_2 <- p21rna[p21rna$variable == cells[example_cell_2_number], c('time', 'value')]

p_p53_cell_1 <- ggplot(p53_example_cell_1, aes(time, value)) +
  geom_line(size = 1.5, color = plot_colors[1]) +
  xlab('Time (h)') + ylab('p53 (a.u.)') +
  custom_plot_theme

p_p21rna_cell_1 <- ggplot(p21rna_example_cell_1, aes(time, value)) +
  geom_line(size = 1.5, color = plot_colors[2]) +
  xlab('Time (h)') + ylab('p21-MS2 (a.u.)') +
  custom_plot_theme

p_p53_cell_2 <- ggplot(p53_example_cell_2, aes(time, value)) +
  geom_line(size = 1.5, color = plot_colors[1]) +
  xlab('Time (h)') + ylab('p53 (a.u.)') +
  custom_plot_theme

p_p21rna_cell_2 <- ggplot(p21rna_example_cell_2, aes(time, value)) +
  geom_line(size = 1.5, color = plot_colors[2]) +
  xlab('Time (h)') + ylab('p21-MS2 (a.u.)') +
  custom_plot_theme

p_example_cells <- 
  plot_grid(
    p_p53_cell_1,
    p_p53_cell_2,
    p_p21rna_cell_1,
    p_p21rna_cell_2,
    ncol = 2,
    labels = c('D', 'F', 'E', 'G'))

################################################################################################################

p53_example_cell_smooth_1 <- smooth_transcription_fractor_dynamics(p53_example_cell_1, plot = FALSE)
p_p53_detrended_1 <- detrend_transcription_factor_dynamics(p53_example_cell_smooth_1, plot = TRUE)
p53_example_cell_smooth_2 <- smooth_transcription_fractor_dynamics(p53_example_cell_2, plot = FALSE)
p_p53_detrended_2 <- detrend_transcription_factor_dynamics(p53_example_cell_smooth_2, plot = TRUE)

p_trends_example_cells <- 
  plot_grid(
    p_p53_detrended_1[[1]],
    p_p53_detrended_2[[1]],
    ncol = 2,
    labels = c('H', 'I'))

################################################################################################################

experiment <- 'IR_2min_resolution'

p53 <- load_data(experiment, 'traces_CFP', FALSE)
p21rna <- load_data(experiment, 'traces_MS2', FALSE)

p_compare_tf_metrics <- compare_predictive_performance_of_metrics_of_transcription(p53, p21rna, detrend_first = TRUE)

################################################################################################################

figure <- 
  plot_grid(
    p_heatmaps,
    plot_grid(
      p_example_cells,
      p_trends_example_cells,
      plot_grid(
        p_compare_tf_metrics,
        NULL,
        ncol = 2),
      ncol = 1,
      labels = c('', '', 'J'),
      rel_heights = c(2, 1, 1)),
    ncol = 1,
    rel_heights = c(0.75, 2))

save_plot('figures/supplemental_figure_p21_transcription_predicted_by_p53_increase.pdf', figure,
          ncol = 2, nrow = 6)
