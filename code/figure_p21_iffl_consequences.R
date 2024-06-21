setwd('/Users/jamiedean/Documents/Papers/Code for Publications/p21-iffl-regulation')

source('code/analysis_functions/determine_modes_of_transcriptional_regulation.R')
source('code/analysis_functions/preprocessing.R')
source('code/analysis_functions/visualize_dynamics.R')

################################################################################################################

p_p21_response_rate_pr <- ggdraw() + 
  draw_image('figures/p21_regulation_ir_pr_dde.png')

p_p21_response_rate_iffl <- ggdraw() + 
  draw_image('figures/p21_regulation_ir_iffl_dde.png')

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

p_p21rna_fano_factor <- 
  ggplot(noise_ratio, aes(time, p21rna_fano_factor, color = condition)) +
    geom_line(size = 2) +
    xlab('Time (h)') + 
    ylab(TeX('$sigma_{p21-MS2}^{2} / mu_{p21-MS2}$')) +
    custom_plot_theme +
    theme(legend.title = element_blank(), legend.justification = c(0, 1), legend.position = c(0, 1),
          legend.box.margin = ggplot2::margin(rep(10, 4)))

p_fano_factor_ratio <-
  ggplot(noise_ratio, aes(time, p21rna_fano_factor/p53_fano_factor, color = condition)) +
    geom_line(size = 2) +
    xlab('Time (h)') + 
    ylab(TeX('$\\frac{sigma_{p21-MS2}^{2} / mu_{p21-MS2}}{sigma_{p53}^{2} / mu_{p53}}$')) +
    custom_plot_theme +
    theme(legend.title = element_blank(), legend.justification = c(0, 1), legend.position = c(0, 1),
          legend.box.margin = ggplot2::margin(rep(10, 4)))

################################################################################################################

# p53 transcriptional targets PR vs IFFL model selection

approach_diagram <- ggdraw() + 
  draw_image('figures/iffl_consequences_diagram.png',
             x = 1, width = 1, height = 1, hjust = 1)

if (file.exists('outputs/rna_seq/rna_seq_transcriptional_regulation_model_simulations')) {
  model_data <- readRDS('outputs/rna_seq/rna_seq_transcriptional_regulation_model_simulations')
} else {
  model_data <- simulate_mrna_dynamics()
}

rf_model <- fit_random_forest_model(model_data)

p53_transcriptional_targets <- as.character(read.csv('data/rna_seq/hafner2017.csv')$gene)

mrna_dynamics_data <- read_rna_seq_data()

mrna_dynamics_p1 <- mrna_dynamics_data$p1[mrna_dynamics_data$p1$gene %in% p53_transcriptional_targets,]
mrna_dynamics_p2 <- mrna_dynamics_data$p2[mrna_dynamics_data$p2$gene %in% p53_transcriptional_targets,]
mrna_dynamics_s1 <- mrna_dynamics_data$s1[mrna_dynamics_data$s1$gene %in% p53_transcriptional_targets,]
mrna_dynamics_s2 <- mrna_dynamics_data$s2[mrna_dynamics_data$s2$gene %in% p53_transcriptional_targets,]

normalized_mrna_dynamics_p1 <- normalize_rna_seq_data(mrna_dynamics_p1, 'fc')
normalized_mrna_dynamics_p2 <- normalize_rna_seq_data(mrna_dynamics_p2, 'fc')
normalized_mrna_dynamics_s1 <- normalize_rna_seq_data(mrna_dynamics_s1, 'fc')
normalized_mrna_dynamics_s2 <- normalize_rna_seq_data(mrna_dynamics_s2, 'fc')

p_model_selection <- model_comparison(rf_model)

transcriptional_regulation_mode_heatmap <- plot_ir_nutlin_versus_ir_mrna_dynamics_heatmap()
png('figures/transcriptional_regulation_mode_heatmap.png', width = 9, height = 10, units = 'in', res = 1200)
transcriptional_regulation_mode_heatmap
dev.off()
p_transcriptional_regulation_mode_heatmap <- ggdraw() + 
  draw_image('figures/transcriptional_regulation_mode_heatmap.png', x = 1, width = 1, height = 1, hjust = 1)

################################################################################################################

# Analysis of TP53 and MDM2 ChIP-seq data from Bevill et al. Cell Genomics 2023

lps141_tp53_mdm2_overlapping_peaks <- find_tp53_mdm2_peak_overlaps('LPS141_untreated')
p_lps141_venn <- ggdraw() + 
  draw_image('figures/LPS141_untreated.png', x = 1, width = 1, height = 1, hjust = 1)
lps853_tp53_mdm2_overlapping_peaks <- find_tp53_mdm2_peak_overlaps('LPS853')
p_lps853_venn <- ggdraw() + 
  draw_image('figures/LPS853.png', x = 1, width = 1, height = 1, hjust = 1)

p_distance_lps141 <- compare_tf_binding_site_to_tss_distances(lps141_tp53_mdm2_overlapping_peaks, 'LPS141')
p_distance_lps853 <- compare_tf_binding_site_to_tss_distances(lps853_tp53_mdm2_overlapping_peaks, 'LPS853')

################################################################################################################

p_microscopy <-
  plot_grid(
    p_p21_response_rate_pr,
    p_p21_response_rate_iffl,
    p_p21rna_fano_factor,
    p_fano_factor_ratio,
    labels = 'AUTO',
    ncol = 2)

p_left <-
  plot_grid(
    p_microscopy,
    approach_diagram,
    labels = c('', 'E'),
    ncol = 1,
    rel_heights = c(1.4, 1)
  )

p_right <-
  plot_grid(
    p_model_selection,
    p_transcriptional_regulation_mode_heatmap,
    labels = c('F', 'G'),
    ncol = 1,
    rel_heights = c(1, 1.75)
  )

p_chip_seq <- 
  plot_grid(
    p_lps141_venn, p_lps853_venn,
    p_distance_lps141, p_distance_lps853,
    labels = c('H', 'I', 'J', 'K'),
    ncol = 4,
    rel_heights = c(1, 1.25))

figure <- 
  plot_grid(
    plot_grid(
      p_left,
      p_right,
      ncol = 2),
    p_chip_seq,
    ncol = 1,
    rel_heights = c(3, 1)
  )

save_plot('figures/figure_p21_iffl_consequences.pdf', figure, ncol = 4.5, nrow = 8)
