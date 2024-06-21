setwd('/Users/jamiedean/Documents/Papers/Code for Publications/p21-iffl-regulation')

source('code/analysis_functions/approximate_bayesian_computation.R')

################################################################################################################

pr_acf_2_5gy <- 
  plot_correlation_function_fit(
    'IR_15min_resolution', 'IR_15min_resolution', '2.5 Gy', 'telegraph_standard', 'p21rna_acf')
pr_ccf_2_5gy <- 
  plot_correlation_function_fit(
    'IR_15min_resolution', 'IR_15min_resolution', '2.5 Gy', 'telegraph_standard', 'p53_p21rna_ccf')
iffl_acf_2_5gy <- 
  plot_correlation_function_fit(
    'IR_15min_resolution', 'IR_15min_resolution', '2.5 Gy', 'telegraph_iffl', 'p21rna_acf')
iffl_ccf_2_5gy <- 
  plot_correlation_function_fit(
    'IR_15min_resolution', 'IR_15min_resolution', '2.5 Gy', 'telegraph_iffl', 'p53_p21rna_ccf')

pr_acf_5gy <- 
  plot_correlation_function_fit(
    'IR_15min_resolution', 'IR_15min_resolution', '5 Gy', 'telegraph_standard', 'p21rna_acf')
pr_ccf_5gy <- 
  plot_correlation_function_fit(
    'IR_15min_resolution', 'IR_15min_resolution', '5 Gy', 'telegraph_standard', 'p53_p21rna_ccf')
iffl_acf_5gy <- 
  plot_correlation_function_fit(
    'IR_15min_resolution', 'IR_15min_resolution', '5 Gy', 'telegraph_iffl', 'p21rna_acf')
iffl_ccf_5gy <- 
  plot_correlation_function_fit(
    'IR_15min_resolution', 'IR_15min_resolution', '5 Gy', 'telegraph_iffl', 'p53_p21rna_ccf')

pr_acf_2min <- 
  plot_correlation_function_fit(
    'IR_2min_resolution', 'IR_2min_resolution', '10 Gy', 'telegraph_standard', 'p21rna_acf')
pr_ccf_2min <- 
  plot_correlation_function_fit(
    'IR_2min_resolution', 'IR_2min_resolution', '10 Gy', 'telegraph_standard', 'p53_p21rna_ccf')
iffl_acf_2min <- 
  plot_correlation_function_fit(
    'IR_2min_resolution', 'IR_2min_resolution', '10 Gy', 'telegraph_iffl', 'p21rna_acf')
iffl_ccf_2min <- 
  plot_correlation_function_fit(
    'IR_2min_resolution', 'IR_2min_resolution', '10 Gy', 'telegraph_iffl', 'p53_p21rna_ccf')

perform_model_selection('IR_15min_resolution', 'IR_15min_resolution', '2.5 Gy',
                        'telegraph_standard', 'telegraph_iffl')
#selected model votes model1 votes model2 post.proba
#             2           48         2952  0.9785167

perform_model_selection('IR_15min_resolution', 'IR_15min_resolution', '5 Gy',
                        'telegraph_standard', 'telegraph_iffl')
#selected model votes model1 votes model2 post.proba
#             2           0         3000  1

perform_model_selection('IR_2min_resolution', 'IR_2min_resolution', '10 Gy',
                        'telegraph_standard', 'telegraph_iffl')
#selected model votes model1 votes model2 post.proba
#             2           0         3000  1

################################################################################################################

figure <- 
  plot_grid(
    pr_acf_2_5gy,
    pr_ccf_2_5gy,
    iffl_acf_2_5gy,
    iffl_ccf_2_5gy,
    pr_acf_5gy,
    pr_ccf_5gy,
    iffl_acf_5gy,
    iffl_ccf_5gy,
    pr_acf_2min,
    pr_ccf_2min,
    iffl_acf_2min,
    iffl_ccf_2min,
    labels = 'AUTO',
    ncol = 4)

save_plot('figures/supplemental_figure_model_comparison.pdf', figure, ncol = 4, nrow = 6)
