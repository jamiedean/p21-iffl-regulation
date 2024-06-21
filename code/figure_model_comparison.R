setwd('/Users/jamiedean/Documents/Papers/Code for Publications/p21-iffl-regulation')

source('code/analysis_functions/approximate_bayesian_computation.R')

################################################################################################################

pr_model_schematic <- ggdraw() + 
  draw_image('figures/positive_regulation_model_schematic.png',
             x = 1, width = 1, height = 1, hjust = 1)

iffl_model_schematic <- ggdraw() + 
  draw_image('figures/incoherent_feedforward_loop_model_schematic.png',
             x = 1, width = 1, height = 1, hjust = 1)

################################################################################################################

experiment_meas <- 'IR_15min_resolution'

p53 <- load_data(experiment_meas, 'traces_CFP', FALSE)
p21rna <- load_data(experiment_meas, 'traces_MS2', FALSE)
cell_number <- 8

params_telegraph_standard <- 
  c(k_on = 0.5,
    k_off = 10000,
    alpha_RNA = 50,
    beta_RNA = 5)

params_telegraph_iffl <- 
  c(k_on = 5,
    k_off = 10000,
    alpha_RNA = 30,
    beta_RNA = 5,
    alpha_Repressor = 1.2,
    beta_Repressor = 1.4)

pr_model_simulation <- 
  plot_single_cell_simulation(p53, p21rna, 'telegraph_standard', params_telegraph_standard, cell_number)
iffl_model_simulation <- 
  plot_single_cell_simulation(p53, p21rna, 'telegraph_iffl', params_telegraph_iffl, cell_number)

model_schematics <- 
  plot_grid(
    pr_model_schematic,
    iffl_model_schematic,
    labels = c('A', 'B'),
    ncol = 1,
    rel_heights = c(1, 1.5))

p_simulations <- 
  plot_grid(
    pr_model_simulation[[1]],
    iffl_model_simulation[[1]],
    labels = c('C', 'D'),
    ncol = 1)

p_autocorrelation_single_cell <- 
  plot_grid(
    pr_model_simulation[[2]],
    iffl_model_simulation[[2]],
    labels = c('E', 'F'),
    ncol = 1)

p_crosscorrelation_single_cell <- 
  plot_grid(
    pr_model_simulation[[3]],
    iffl_model_simulation[[3]],
    labels = c('G', 'H'),
    ncol = 1)

p_top_row <- 
  plot_grid(
    model_schematics,
    p_simulations,
    p_autocorrelation_single_cell,
    p_crosscorrelation_single_cell,
    ncol = 4)

################################################################################################################

experiment_meas <- 'IR_15min_resolution'
experiment_pred <- 'IR_15min_resolution'

pr_acf_10gy <- 
  plot_correlation_function_fit(experiment_meas, experiment_pred, '10 Gy', 'telegraph_standard', 'p21rna_acf')
pr_ccf_10gy <- 
  plot_correlation_function_fit(experiment_meas, experiment_pred, '10 Gy', 'telegraph_standard', 'p53_p21rna_ccf')
iffl_acf_10gy <- 
  plot_correlation_function_fit(experiment_meas, experiment_pred, '10 Gy', 'telegraph_iffl', 'p21rna_acf')
iffl_ccf_10gy <- 
  plot_correlation_function_fit(experiment_meas, experiment_pred, '10 Gy', 'telegraph_iffl', 'p53_p21rna_ccf')

p_correlation_functions <-
  plot_grid(
    pr_acf_10gy,
    pr_ccf_10gy,
    iffl_acf_10gy,
    iffl_ccf_10gy,
    labels = c('I', 'J', 'K', 'L'),
    ncol = 4)

perform_model_selection(experiment_meas, experiment_pred, '10 Gy', 'telegraph_standard', 'telegraph_iffl')
#selected model votes model1 votes model2 post.proba
#             2           0         3000  1

################################################################################################################

figure <- 
  plot_grid(
    p_top_row,
    p_correlation_functions,
    ncol = 1,
    rel_heights = c(1, 0.5))

save_plot('figures/figure_model_comparison.pdf', figure, ncol = 4, nrow = 4)
