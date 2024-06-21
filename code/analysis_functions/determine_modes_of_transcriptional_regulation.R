library(circlize)
library(ComplexHeatmap)
library(deSolve)
library(MatrixGenerics)
library(randomForest)
library(splines)

source('code/analysis_functions/chip_seq_analysis.R')
source('code/analysis_functions/custom_plot_theme.R')
source('code/analysis_functions/rna_seq_data_preprocessing.R')

################################################################################################################

smooth_p53_dynamics <- function(plot = FALSE) {
  
  # p53 protein expression data from Hafner et al. 2017 Nat Struct Mol Biol
  normalized_p53_dynamics_p <- 
    data.frame(time = seq(0, 12), 
               p53 = c(1, 4.37, 12.78, 5.20, 1.98, 1.50, 1.11, 3.20, 4.01, 1.28, 1.43, 1.85, 2.01))
  normalized_p53_dynamics_s <-
    data.frame(time = seq(0, 12), 
               p53 = c(1, 4.37, 12.78, 15.74, 25.69, 22.41, 21.50, 36.88, 25.02, 37.43, 29.85, 15.41, 22.63))
  
  p53_p_spline <- smooth.spline(normalized_p53_dynamics_p, nknots = 12)
  p53_p <- predict(p53_p_spline, seq(0, 12, 0.1))
  p53_s_spline <- smooth.spline(normalized_p53_dynamics_s, nknots = 12)
  p53_s <- predict(p53_s_spline, seq(0, 12, 0.1))
  
  if (plot == TRUE) {
    plot(normalized_p53_dynamics_p, type = 'l', ylim = c(0, 50))
    lines(normalized_p53_dynamics_s)
    lines(p53_p, col = plot_colors[1])
    lines(p53_s, col = plot_colors[2]) 
  }
  
  p53_dynamcis_smooth <- data.frame(p53_p = p53_p, p53_s = p53_s)
  
  return(p53_dynamcis_smooth)
}

mrna_dynamics_ode_model <-function(t, x, params, p53_dynamics) {
  
  p53_p_t <- p53_dynamics$p53_p.y[which.min(abs(p53_dynamics$p53_p.x - t))]
  p53_s_t <- p53_dynamics$p53_s.y[which.min(abs(p53_dynamics$p53_s.x - t))]
  mrna_p <- x[1]
  mrna_s <- x[2]
  mrna_iffl_s <- x[3]
  
  alpha <- params['alpha']
  beta <- params['beta']
  K <- params['K']
  gamma <- params['gamma']

  gamma_t <- ifelse(t > 3, gamma, 0)
  
  dmrna_p_dt <- alpha*p53_p_t/(K + p53_p_t) - beta*mrna_p
  dmrna_s_dt <- alpha*p53_s_t/(K + p53_s_t) - beta*mrna_s
  dmrna_s_iffl_dt <- (alpha + gamma_t)*p53_s_t/(K + p53_s_t) - beta*mrna_s
  
  list(c(dmrna_p_dt, dmrna_s_dt, dmrna_s_iffl_dt))
}

simulate_mrna_dynamics <- function(n_genes = 10000) {
  
  params <- data.frame(alpha = rlnorm(n_genes/2, 1, 1.5),
                       K = abs(rnorm(n_genes/2, 0.1, 10)),
                       beta = abs(rnorm(n_genes/2, 0, 0.3)),
                       gamma = abs(rnorm(n_genes/2, 1, 10)))
  
  times <- seq(0, 12)
  
  p53_dynamics_smooth <- smooth_p53_dynamics()
  
  sim_mrna_p <- matrix(nrow = dim(params)[1], ncol = length(times))
  sim_mrna_s <- matrix(nrow = dim(params)[1], ncol = length(times))
  sim_mrna_s_iffl <- matrix(nrow = dim(params)[1], ncol = length(times))
  
  for (i in 1:dim(params)[1]) {
    
    y_init <- c(mrna_p = 1, mrna_s = 1, mrna_s_iffl = 1)
    out <- ode(func = mrna_dynamics_ode_model, y = y_init, times = times, parms = as.matrix(params)[i,],
               p53_dynamics = p53_dynamics_smooth)
    sim_mrna_p[i,] <- out[, 'mrna_p']
    sim_mrna_s[i,] <- out[, 'mrna_s']
    sim_mrna_s_iffl[i,] <- out[, 'mrna_s_iffl']
  }
  
  sim_mrna_data_p <- rbind(sim_mrna_p, sim_mrna_p)
  sim_mrna_data_s <- rbind(sim_mrna_s, sim_mrna_s_iffl)
  
  heatmap_data <- cbind(sim_mrna_data_p, sim_mrna_data_s)
  
  gamma <- c(rep(0, dim(params)[1]), params$gamma)
  
  model_data <- data.frame(gamma, sum_mrna_p = rowSums(sim_mrna_data_p), sum_mrna_s = rowSums(sim_mrna_data_s),
                           sum_mrna_ratio = rowSums(sim_mrna_data_s/sim_mrna_data_p))
  
  model_data <- data.frame(cbind(
    sim_mrna_data_p, sim_mrna_data_s,
    c(rep(1, dim(sim_mrna_data_p)[1]/2), rep(2, dim(sim_mrna_data_p)[1]/2))))
  colnames(model_data) <- c(paste0('p', seq(0, 12)), paste0('s', seq(0, 12)), 'model_number')
  
  saveRDS(model_data, 'outputs/rna_seq/rna_seq_transcriptional_regulation_model_simulations')
  
  return(model_data)
}

fit_random_forest_model <- function(model_data) {
  
  rf_model <- randomForest(as.factor(model_number) ~ ., data = model_data, proximity = TRUE)
  print(rf_model)
  round(importance(rf_model), 2)
  
  return(rf_model)
}

model_comparison <- function(rf_model) {
  
  normalized_mrna_dynamics_p_mean <- 
    apply(simplify2array(list(normalized_mrna_dynamics_p1, normalized_mrna_dynamics_p2)), 1:2, mean)
  normalized_mrna_dynamics_s_mean <- 
    apply(simplify2array(list(normalized_mrna_dynamics_s1, normalized_mrna_dynamics_s2)), 1:2, mean)
  
  pred_data <- data.frame(cbind(normalized_mrna_dynamics_p_mean, normalized_mrna_dynamics_s_mean))
  colnames(pred_data) <- c(paste0('p', seq(0, 12)), paste0('s', seq(0, 12)))
  rownames(pred_data) <- rownames(normalized_mrna_dynamics_p_mean)
  
  pred_data$model_number <- c(t(predict(rf_model, pred_data)))
  pred_data$model_number
  pred_data$model_name <- ifelse(pred_data$model_number == 2, 'IFFL Model', 'PR Model')
  pred_data$prob <- c(t(predict(rf_model, pred_data, type = 'prob')[, 2]))
  pred_data$gene <- rownames(pred_data)
  pred_data[which(pred_data$model_name == 'IFFL Model'),]
  hist(pred_data$prob, breaks = 30)
  pred_data[which(pred_data$prob > 0.75),]
  
  model_data$model_name <- ifelse(model_data$model_number == 2, 'IFFL Model', 'PR Model')
  model_data$sum_mrna_p <- rowSums(model_data[, 1:13])
  model_data$sum_mrna_s <- rowSums(model_data[, 14:26])
  model_data$sum_mrna_ratio <- rowSums(model_data[, 14:26]/model_data[, 1:13])
  model_data$model_name <- factor(model_data$model_name, levels = c('PR Model', 'IFFL Model'))
  
  pred_data$sum_mrna_p <- rowSums(normalized_mrna_dynamics_p_mean)
  pred_data$sum_mrna_s <- rowSums(normalized_mrna_dynamics_s_mean)
  pred_data$sum_mrna_ratio <- rowSums(normalized_mrna_dynamics_s_mean/normalized_mrna_dynamics_p_mean)
  pred_data[pred_data$prob < 0.75, 'gene'] <- ''
  
  p <- ggplot() +
    stat_density2d(aes(log(sum_mrna_p), log(sum_mrna_ratio), fill = model_name, alpha = ..level..),
                   data = model_data, geom = 'polygon') +
    geom_point(aes(log(sum_mrna_p), log(sum_mrna_ratio), col = 'Data'),
               data = pred_data) +
    geom_text_repel(aes(log(sum_mrna_p), log(sum_mrna_ratio), label = gene),
                    data = pred_data, size = 5, box.padding = 1) +
    scale_fill_manual(values = c('PR Model' = plot_colors[1], 'IFFL Model' = plot_colors[2])) +
    scale_color_manual(values = c('Data' = plot_colors[3])) +
    xlab(TeX('log $\\left(\\sum mRNA_{IR}(t)/mRNA_{IR}(0) \\right)$')) + 
    ylab(TeX(
      'log $\\left(\\sum \\frac{mRNA_{IR+Nut}(t)/mRNA_{IR+Nut}(0)}{mRNA_{IR}(t)/mRNA_{IR}(0)} \\right)$')) + 
    theme(legend.position = c(0.7, 0.75),
          legend.direction = 'vertical',
          legend.box = 'horizontal',
          legend.spacing.x = unit(0.1, 'cm'))
  
  print(pred_data)
  write.csv(pred_data$prob, 'outputs/rna_seq/predicted_transcription_mode_p53_targets.csv',
            row.names = row.names(pred_data),
            col.names = c('gene', 'prob_iffl'))
  
  return(p)
}

plot_ir_nutlin_versus_ir_mrna_dynamics_heatmap <- function() {
  
  normalized_mrna_dynamics_diff_mean <-
    apply(simplify2array(list(normalized_mrna_dynamics_s1/normalized_mrna_dynamics_p1,
                              normalized_mrna_dynamics_s2/normalized_mrna_dynamics_p2)), 1:2, mean)
  
  heatmap_data <- normalized_mrna_dynamics_diff_mean
  colnames(heatmap_data) <- as.character(seq(0, 12))
  
  log_fc_mrna_dynamics_p1 <- normalize_rna_seq_data(mrna_dynamics_p1, 'log_fc')
  log_fc_mrna_dynamics_p2 <- normalize_rna_seq_data(mrna_dynamics_p2, 'log_fc')
  log_fc_mrna_dynamics_s1 <- normalize_rna_seq_data(mrna_dynamics_s1, 'log_fc')
  log_fc_mrna_dynamics_s2 <- normalize_rna_seq_data(mrna_dynamics_s2, 'log_fc')
  
  log_fc_mrna_dynamics_p_mean <- 
    apply(simplify2array(list(log_fc_mrna_dynamics_p1, log_fc_mrna_dynamics_p2)), 1:2, mean)
  log_fc_mrna_dynamics_s_mean <- 
    apply(simplify2array(list(log_fc_mrna_dynamics_s1, log_fc_mrna_dynamics_s2)), 1:2, mean)
  
  log2_basal <- log2(c(t(colMeans(rbind(mrna_dynamics_p1[, 1], mrna_dynamics_p2[, 1])))))
  fold_ir <- rowMaxs(log_fc_mrna_dynamics_p_mean)
  fold_ir_nutlin <- rowMaxs(log_fc_mrna_dynamics_s_mean)
  
  p53_chip_p53_targets_max_peak <- select_p53_targets_chip_seq_max_peaks()
  p53_chip_p53_targets_max_peak_basal <- 
    log2(p53_chip_p53_targets_max_peak$NormReads_p53_t0/p53_chip_p53_targets_max_peak$NormReads_input_t0)
  p53_chip_p53_targets_max_peak_tss_distance <- log(abs(p53_chip_p53_targets_max_peak$DisttoTSS))
  p53_chip_p53_targets_max_peak_hmm <- p53_chip_p53_targets_max_peak$chromatin_hmm
  p53_chip_p53_targets_max_peak_hmm_tss <- p53_chip_p53_targets_max_peak$chromatin_hmm_tss
  
  prob_iffl <- read.csv('outputs/rna_seq/predicted_transcription_mode_p53_targets.csv')$x
  
  row_ann <- 
    HeatmapAnnotation(
      'P(IFFL)' = prob_iffl,
      'IFFL' = prob_iffl > 0.5,
      'log TSS Distance' = p53_chip_p53_targets_max_peak_tss_distance,
      which = 'row',
      annotation_name_side = 'top',
      col = list(
        'P(IFFL)' = colorRamp2(quantile(prob_iffl, c(0, 1)), c('white', plot_colors[2])),
        'IFFL' = structure(c('white', 'black'), names = c('FALSE', 'TRUE')),
        'log TSS Distance' = 
          colorRamp2(quantile(p53_chip_p53_targets_max_peak_tss_distance, c(0, 1)), c('white', 'pink'))),
      border = TRUE,
      annotation_name_gp = gpar(fontsize = 14),
      annotation_legend_param = list(
        title_gp = gpar(fontsize = 14),
        labels_gp = gpar(fontsize = 14)))
  
  ht <- Heatmap(
    heatmap_data,
    name = 'IR + Nut/IR',
    col = colorRamp2(quantile(heatmap_data, c(0, 1)), c('white', plot_colors[1])),
    clustering_method_rows = 'ward.D2',
    row_dend_reorder = TRUE,
    row_dend_side = 'left',
    row_split = 3,
    cluster_columns = FALSE,
    row_names_side = 'right',
    column_names_side = 'top',
    column_names_rot = 0,
    column_title = 'Time (h)',
    border = TRUE,
    right_annotation = row_ann,
    column_title_gp = gpar(fontsize = 14),
    row_title_gp = gpar(fontsize = 14),
    column_names_gp = gpar(fontsize = 14),
    row_names_gp = gpar(fontsize = 14),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 14),
      labels_gp = gpar(fontsize = 14)))
  
  return(ht)
}
