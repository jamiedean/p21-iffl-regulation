setwd('/Users/jamiedean/Documents/Papers/Code for Publications/p21-iffl-regulation')

source('code/analysis_functions/prediction_of_transcription_by_transcription_factor_metrics.R')

################################################################################################################

experiment <- 'IR_15min_resolution'

p53 <- load_data(experiment, 'traces_CFP', FALSE)
p21rna <- load_data(experiment, 'traces_MS2', FALSE)
p21 <- load_data(experiment, 'traces_Texas', FALSE)

cells <- intersect(unique(p53$variable), unique(p21rna$variable))

example_cell_number <- 40
p53_example_cell <- p53[p53$variable == cells[example_cell_number], c('time', 'value')]
p21rna_example_cell <- p21rna[p21rna$variable == cells[example_cell_number], c('time', 'value')]
p21_example_cell <- p21[p21$variable == cells[example_cell_number], c('time', 'value')]

p_p53 <- ggplot(p53_example_cell, aes(time, value)) +
  geom_line(linewidth = 1.5, color = plot_colors[4]) +
  xlab('Time (h)') + ylab('p53 (a.u.)') +
  custom_plot_theme

p_p21rna <- ggplot(p21rna_example_cell, aes(time, value)) +
  geom_line(linewidth = 1.5, color = plot_colors[4]) +
  xlab('Time (h)') + ylab('p21-MS2 (a.u.)') +
  custom_plot_theme

p_p21 <- ggplot(p21_example_cell, aes(time, value)) +
  geom_line(linewidth = 1.5, color = plot_colors[4]) +
  xlab('Time (h)') + ylab('p21 (a.u.)') +
  custom_plot_theme

p_p53_level <- smooth_transcription_fractor_dynamics(p53_example_cell, plot = TRUE)

p53_example_cell_smooth <- smooth_transcription_fractor_dynamics(p53_example_cell, plot = FALSE)
p_p53_detrended <- detrend_transcription_factor_dynamics(p53_example_cell_smooth, plot = TRUE, num_metrics = 2)
p_p53_change <- calculate_change_in_tf_level(p53_example_cell_smooth, plot = TRUE, increase_only = FALSE)

p_gene_state_example_cell <- infer_gene_state(p21rna_example_cell, plot = TRUE, num_metrics = 2)

p_logistic_regression_p53_level <- 
  predict_promoter_status(p53, p21rna, example_cell_number, 'level', num_metrics = 2, TRUE, TRUE)
p_logistic_regression_p53_change <- 
  predict_promoter_status(p53, p21rna, example_cell_number, 'change', num_metrics = 2, TRUE, TRUE)

p_compare_tf_metrics_auc <- 
  compare_predictive_performance_of_metrics_of_transcription(p53, p21rna, num_metrics = 2, detrend_first = TRUE)

figure <- plot_grid(
  p_p53,
  p_p21rna,
  p_p21,
  p_p53_detrended[[2]],
  p_p53_change,
  p_gene_state_example_cell,
  p_logistic_regression_p53_level,
  p_logistic_regression_p53_change,
  p_compare_tf_metrics_auc,
  labels = 'AUTO',
  ncol = 3
)

save_plot('figures/figure_p21_transcription_predicted_by_p53_increase.pdf', figure, ncol = 2, nrow = 3)
