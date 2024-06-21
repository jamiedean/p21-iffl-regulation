setwd('/Users/jamiedean/Documents/Papers/Code for Publications/p21-iffl-regulation')

source('code/analysis_functions/approximate_bayesian_computation.R')
source('code/analysis_functions/prediction_of_transcription_by_transcription_factor_metrics.R')

################################################################################################################

experiment_meas <- 'IR_Nutlin'
experiment_pred <- 'IR_Nutlin'

pr_acf_ir <- 
  plot_correlation_function_fit(experiment_meas, experiment_pred, 'IR', 'telegraph_standard', 'p21rna_acf')
pr_ccf_ir <- 
  plot_correlation_function_fit(experiment_meas, experiment_pred, 'IR', 'telegraph_standard', 'p53_p21rna_ccf')
iffl_acf_ir <- 
  plot_correlation_function_fit(experiment_meas, experiment_pred, 'IR', 'telegraph_iffl', 'p21rna_acf')
iffl_ccf_ir <- 
  plot_correlation_function_fit(experiment_meas, experiment_pred, 'IR', 'telegraph_iffl', 'p53_p21rna_ccf')

pr_acf_ir_nutlin <- 
  plot_correlation_function_fit(experiment_meas, experiment_pred, 'IR + nutlin', 'telegraph_standard', 'p21rna_acf')
pr_ccf_ir_nutlin <- 
  plot_correlation_function_fit(experiment_meas, experiment_pred, 'IR + nutlin', 'telegraph_standard', 'p53_p21rna_ccf')
iffl_acf_ir_nutlin <- 
  plot_correlation_function_fit(experiment_meas, experiment_pred, 'IR + nutlin', 'telegraph_iffl', 'p21rna_acf')
iffl_ccf_ir_nutlin <- 
  plot_correlation_function_fit(experiment_meas, experiment_pred, 'IR + nutlin', 'telegraph_iffl', 'p53_p21rna_ccf')

perform_model_selection(experiment_meas, experiment_pred, 'IR', 'telegraph_standard', 'telegraph_iffl')
#selected model votes model1 votes model2 post.proba
#             2            0         3000          1

perform_model_selection(experiment_meas, experiment_pred, 'IR + nutlin', 'telegraph_standard', 'telegraph_iffl')
#selected model votes model1 votes model2 post.proba
#             2         1030         1970  0.6978167

################################################################################################################

experiment <- 'IR_Nutlin'

p53 <- load_data(experiment, 'traces_CFP', FALSE)
p21rna <- load_data(experiment, 'traces_YFP', FALSE)

p53_ir <- p53[p53$condition == 'IR',]
p21rna_ir <- p21rna[p21rna$condition == 'IR',]
p53_ir_nutlin <- p53[p53$condition == 'IR + nutlin',]
p21rna_ir_nutlin <- p21rna[p21rna$condition == 'IR + nutlin',]

p_compare_tf_metrics_ir <- 
  compare_predictive_performance_of_metrics_of_transcription(
    p53_ir, p21rna_ir, num_metrics = 2, detrend_first = TRUE)
p_compare_tf_metrics_ir_nutlin <- 
  compare_predictive_performance_of_metrics_of_transcription(
    p53_ir_nutlin, p21rna_ir_nutlin, num_metrics = 2, detrend_first = FALSE)

################################################################################################################

title_ir <- ggdraw() +
  draw_label('IR', fontface = 'bold', size = 22)

title_ir_nutlin <- ggdraw() +
  draw_label('IR + nutlin', fontface = 'bold', size = 22)

plot_ir <- plot_grid(
  title_ir,
  plot_grid(
    plot_grid(
      pr_acf_ir,
      pr_ccf_ir,
      iffl_acf_ir,
      iffl_ccf_ir,
      labels = c('A', 'B', 'E', 'F'),
      ncol = 2),
    plot_grid(
      NULL, p_compare_tf_metrics_ir, NULL,
      labels = c('', 'I', ''),
      ncol = 3,
      rel_widths = c(0.25, 0.5, 0.25)
    ),
    ncol = 1,
    rel_heights = c(2, 1)),
  ncol = 1,
  rel_heights = c(1, 10))

plot_ir_nutlin <- plot_grid(
  title_ir_nutlin,
  plot_grid(
    plot_grid(
      pr_acf_ir_nutlin,
      pr_ccf_ir_nutlin,
      iffl_acf_ir_nutlin,
      iffl_ccf_ir_nutlin,
      labels = c('C', 'D', 'G', 'H'),
      ncol = 2),
    plot_grid(
      NULL, p_compare_tf_metrics_ir_nutlin, NULL,
      labels = c('', 'J', ''),
      ncol = 3,
      rel_widths = c(0.25, 0.5, 0.25)
    ),
    ncol = 1,
    rel_heights = c(2, 1)),
  ncol = 1,
  rel_heights = c(1, 10))

figure <- plot_grid(
  plot_ir,
  plot_ir_nutlin,
  ncol = 2)

save_plot('figures/figure_model_validation.pdf', figure, ncol = 3.5, nrow = 4.5)
