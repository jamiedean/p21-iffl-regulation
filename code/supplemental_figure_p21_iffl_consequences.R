setwd('/Users/jamiedean/Documents/Papers/Code for Publications/p21-iffl-regulation')

source('code/analysis_functions/determine_modes_of_transcriptional_regulation.R')
source('code/analysis_functions/preprocessing.R')
source('code/analysis_functions/visualize_dynamics.R')

################################################################################################################

experiment <- 'IR_Nutlin'

p53 <- load_data(experiment, 'traces_CFP', FALSE)
p21rna <- load_data(experiment, 'traces_YFP', FALSE)

p53_dynamics <- plot_summary_dynamics_by_condition(p53, experiment, log_transform_expression = TRUE)
p21rna_dynamics <- plot_summary_dynamics_by_condition(p21rna, experiment)

noise_ratio <- 
  data.frame(time = p53_dynamics$time,
             p53_mean = p53_dynamics$mean,
             p21rna_mean = p21rna_dynamics$mean,
             p53_coef_var = p53_dynamics$coef_var,
             p21rna_coef_var = p21rna_dynamics$coef_var,
             p53_fano_factor = p53_dynamics$fano_factor,
             p21rna_fano_factor = p21rna_dynamics$fano_factor,
             condition = p53_dynamics$condition)

p_p21rna_coef_var <-
  ggplot(noise_ratio, aes(time, p21rna_coef_var, color = condition)) +
  geom_line(size = 2) +
  xlab('Time (h)') + 
  ylab(TeX('$sigma_{p21-MS2} / mu_{p21-MS2}$')) +
  custom_plot_theme +
  theme(legend.title = element_blank(), legend.justification = c(1, 1), legend.position = c(1, 1),
        legend.box.margin = ggplot2::margin(rep(10, 4)))

p_p53_fano_factor <-
  ggplot(noise_ratio, aes(time, p53_fano_factor, color = condition)) +
  geom_line(size = 2) +
  xlab('Time (h)') + 
  ylab(TeX('$sigma_{p53}^{2} / mu_{p53}$')) +
  custom_plot_theme +
  theme(legend.title = element_blank(), legend.justification = c(1, 1), legend.position = c(1, 1),
        legend.box.margin = ggplot2::margin(rep(10, 4)))

p_coef_var_ratio <-
  ggplot(noise_ratio, aes(time, p21rna_coef_var/p53_coef_var, color = condition)) +
  geom_line(size = 2) +
  xlab('Time (h)') + 
  ylab(TeX('$\\frac{sigma_{p21-MS2} / mu_{p21-MS2}}{sigma_{p53} / mu_{p53}}$')) +
  custom_plot_theme +
  theme(legend.title = element_blank(), legend.justification = c(1, 1), legend.position = c(1, 1),
        legend.box.margin = ggplot2::margin(rep(10, 4)))

################################################################################################################

experiment <- 'IR_15min_resolution'

p53 <- load_data(experiment, 'traces_CFP', FALSE)
p21rna <- load_data(experiment, 'traces_MS2', FALSE)

p53_dynamics <- plot_summary_dynamics_by_condition(p53, experiment, log_transform_expression = TRUE)
p21rna_dynamics <- plot_summary_dynamics_by_condition(p21rna, experiment)

noise_ratio <- 
  data.frame(time = p53_dynamics$time,
             p53_mean = p53_dynamics$mean,
             p21rna_mean = p21rna_dynamics$mean,
             p53_coef_var = p53_dynamics$coef_var,
             p21rna_coef_var = p21rna_dynamics$coef_var,
             p53_fano_factor = p53_dynamics$fano_factor,
             p21rna_fano_factor = p21rna_dynamics$fano_factor,
             condition = p53_dynamics$condition)

p_p21rna_fano_factor <- 
  ggplot(noise_ratio, aes(time, p21rna_fano_factor, color = condition)) +
    geom_line(size = 2) +
    xlab('Time (h)') + 
    ylab(TeX('$sigma_{p21-MS2}^{2} / mu_{p21-MS2}$')) +
    custom_plot_theme +
    theme(legend.title = element_blank(), legend.justification = c(1, 1), legend.position = c(1, 1),
          legend.box.margin = ggplot2::margin(rep(10, 4)))

p_fano_factor_ratio <-
  ggplot(noise_ratio, aes(time, p21rna_fano_factor/p53_fano_factor, color = condition)) +
    geom_line(size = 2) +
    xlab('Time (h)') + 
    ylab(TeX('$\\frac{sigma_{p21-MS2}^{2} / mu_{p21-MS2}}{sigma_{p53}^{2} / mu_{p53}}$')) +
    custom_plot_theme +
    theme(legend.title = element_blank(), legend.justification = c(1, 1), legend.position = c(1, 1),
          legend.box.margin = ggplot2::margin(rep(10, 4)))

################################################################################################################

# Analysis of TP53 and MDM2 ChIP-seq data from Bevill et al. Cell Genomics 2023

hct116_tp53_mdm2_overlapping_peaks <- find_tp53_mdm2_peak_overlaps('HCT116')
p_hct116_venn <- ggdraw() + 
  draw_image('figures/HCT116.png', x = 1, width = 1, height = 1, hjust = 1)
u2os_tp53_mdm2_overlapping_peaks <- find_tp53_mdm2_peak_overlaps('U2OS')
p_u2os_venn <- ggdraw() + 
  draw_image('figures/U2OS.png', x = 1, width = 1, height = 1, hjust = 1)

################################################################################################################

figure <- 
  plot_grid(
    p_p21rna_coef_var,
    p_p53_fano_factor,
    p_coef_var_ratio,
    p_p21rna_fano_factor,
    p_fano_factor_ratio,
    NULL,
    p_hct116_venn,
    p_u2os_venn,
    labels = c('A', 'B', 'C', 'D', 'E', '', 'F', 'G'),
    ncol = 3)

save_plot('figures/supplemental_figure_p21_iffl_consequences.pdf', figure, ncol = 3, nrow = 4)
