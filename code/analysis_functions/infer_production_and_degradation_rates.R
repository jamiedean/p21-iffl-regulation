library(dpseg)
library(ggstance)
library(mclust)
library(nlraa)

source('code/analysis_functions/custom_plot_theme.R')
source('code/analysis_functions/preprocessing.R')

################################################################################################################

infer_protein_degradation_rates_single_cell <- function(traces, cell_number, start_frame = 110) {
  
  cell <- unique(traces$variable)[cell_number]
  cell_time <- traces[traces$variable == cell, 'time'][start_frame:length(unique(traces$time))]
  cell_fluorescence <- traces[traces$variable == cell, 'value'][start_frame:length(unique(traces$time))]
  
  piecewise_regression <- dpseg(x = cell_time, y = log(cell_fluorescence), P = 0.01, store.matrix = TRUE)
  
  regression_line_1 <- 
    data.frame(
      x = seq(piecewise_regression$segments$x1[1] - 0.2, piecewise_regression$segments$x2[1] + 0.2, 0.1),
      y = piecewise_regression$segments$intercept[1] + piecewise_regression$segments$slope[1]*
        seq(piecewise_regression$segments$x1[1] - 0.2, piecewise_regression$segments$x2[1] + 0.2, 0.1))
  
  regression_line_2 <- 
    data.frame(
      x = seq(piecewise_regression$segments$x1[2] - 0.2, piecewise_regression$segments$x2[2] +0.2, 0.1),
      y = piecewise_regression$segments$intercept[2] + piecewise_regression$segments$slope[2]*
        seq(piecewise_regression$segments$x1[2] - 0.2, piecewise_regression$segments$x2[2] + 0.2, 0.1))
  
  p_infer_degradation_rate <-
    ggplot() +
    geom_line(data = data.frame(cell_time, cell_fluorescence), aes(cell_time, log(cell_fluorescence)),
              color = plot_colors[4], linewidth = 2) +
    geom_line(data = regression_line_1, aes(x, y), color = plot_colors[1], linewidth = 2, alpha = 0.8) +
    geom_line(data = regression_line_2, aes(x, y), color = plot_colors[2], linewidth = 2, alpha = 0.8) +
    xlab('Time (h)') + ylab('log(p21) (a.u.)') +
    custom_plot_theme
  
  return(p_infer_degradation_rate)
}

infer_protein_degradation_rates <- function(traces, start_frame = 110) {
  # Note that as p21 production is inhibited through addition of p53 and p21 siRNA cells transition
  # through G1/S leading to rapid degradation of p21 resulting in a 2 phase decay (slow followed by fast)
  # Therefore, degradation rate is inferred through piecewise linear regression (on the log scale)
  
  # Cells were exposed to 40nM siRNA targeting p53 alone, p53 + p21 or control at frame 80 (20 h)
  
  degradation_rates <- c()
  
  for (cell_number in 1:length(unique(traces$variable))) {
    
    cell <- unique(traces$variable)[cell_number]
    cell_time <- traces[traces$variable == cell, 'time'][start_frame:length(unique(traces$time))]
    cell_fluorescence <- traces[traces$variable == cell, 'value'][start_frame:length(unique(traces$time))]
    
    piecewise_regression <- dpseg(x = cell_time, y = log(cell_fluorescence), P = 0.0001, store.matrix = TRUE)

    degradation_rates <- append(degradation_rates, max(-piecewise_regression$segments$slope))
  }
  
  degradation_rates <- degradation_rates[degradation_rates > 0]
  
  mixture_model <- Mclust(log(degradation_rates), G = 2, model = 'V', na.rm = TRUE)
  mixture_model_density <- densityMclust(log(degradation_rates), plot = FALSE)
  
  inferred_degradation_rates <- exp(mixture_model$parameters$mean)
  
  degradation_rate_df <- 
    data.frame(log_degradation_rate = log(degradation_rates), mixture_model = mixture_model_density$density)
  
  p_degradation_rates <- ggplot(degradation_rate_df) +
    geom_histogram(aes(x = log_degradation_rate, y = after_stat(density)),
                   bins = 50, fill = plot_colors[4], color = 'white', alpha = 0.5) +
    geom_line(aes(x = log_degradation_rate, y = mixture_model), colour = plot_colors[4], linewidth = 1.5) +
    geom_vline(xintercept = log(inferred_degradation_rates[1]), linewidth = 1.5, colour = plot_colors[1]) +
    geom_text(aes(x = log(inferred_degradation_rates[1]) - 0.9, 
                  label = paste(round(inferred_degradation_rates[1], 2), '~h^{-1}'),
                  y = 0.65),
              color = plot_colors[1], size = 6, parse = TRUE) +
    geom_vline(xintercept = log(inferred_degradation_rates[2]), linewidth = 1.5, colour = plot_colors[2]) +
    geom_text(aes(x = log(inferred_degradation_rates[2]) - 0.9,
                  label = paste(round(inferred_degradation_rates[2], 2), '~h^{-1}'),
                  y = 0.65),
              color = plot_colors[2], size = 6, parse = TRUE) +
    ylim(c(0, 0.7)) + 
    xlab(TeX('log(p21 degradation rate) (h$^{-1}$)')) + ylab('Density') +
    custom_plot_theme
  
  return(list(p_degradation_rates, mixture_model$parameters))
}

determine_steady_state_cells <- function(traces, time_from_end_of_experiment) {
  # Fit linear regression to end of protein time course data for individual cells to determine if protein expression has
  # reached steady state
  
  start_time <- max(traces$time) - time_from_end_of_experiment
  end_time <- max(traces$time)
  
  linear_regression_fits <- data.frame(matrix(nrow = length(unique(traces$variable)), ncol = 2))
  colnames(linear_regression_fits) <- c('slope', 'p_value')
  rownames(linear_regression_fits) <- unique(traces$variable)
  
  for (cell_number in 1:length(unique(traces$variable))) {
    
    cell <- unique(traces$variable)[cell_number]
    traces_cell <- traces[traces$variable == cell,] 
    regression_model <- lm(log(value) ~ time,
                           traces_cell[traces_cell$time >= start_time & traces_cell$time <= end_time,])
    linear_regression_fits[cell_number, 'slope'] <- summary(regression_model)$coefficients['time', 'Estimate']
    linear_regression_fits[cell_number, 'p_value'] <- summary(regression_model)$coefficients['time', 'Pr(>|t|)']
  }
  
  p_determine_steady_state <- 
    ggplot(linear_regression_fits, aes(slope, -log10(p_value), color = factor(abs(slope) < 0.1))) +
    geom_point() +
    custom_plot_theme
  print(p_determine_steady_state)
  
  steady_state_cells <- rownames(linear_regression_fits[which(abs(linear_regression_fits$slope) < 0.05),])
  
  return(steady_state_cells)
}

infer_production_rates <- function(experiment, p21_degradation_rate_mixture_model_fit, time_from_end_of_experiment = 2.5) {
  # dp21/dt = kp*p53/(K + p53) - kd*p21
  # dp21/dt + kd*p21 = kp*p53/(K + p53)
  # If dp21/dt = 0 then,
  # p21/p53 = kp/(kd*(K + p53))
  # log(p21/p53) = log(kp/kd) - log(K + p53)
  # If p53 >> K then,
  # log(p21/p53) ~ log(kp/kd) - log(p53)
  # log(p21) - log(p53) ~ log(kp/kd) - log(p53)
  # log(p21) ~ log(kp/kd)
  
  fast_degradation_sigma <- sqrt(p21_degradation_rate_mixture_model_fit$variance$sigmasq[1])
  slow_degradation_sigma <- sqrt(p21_degradation_rate_mixture_model_fit$variance$sigmasq[2])
  
  p53 <- load_data(experiment, 'traces_CFP', FALSE)
  p21 <- load_data(experiment, 'traces_Texas', FALSE)
  
  # Only include cells where p21 protein expression is approximately at steady state 
  p21_steady_state_cells <- determine_steady_state_cells(p21, time_from_end_of_experiment)
  p21 <- p21[p21$variable %in% p21_steady_state_cells,]
  
  # Use only data points from the end of protein time course so that protein expression is approximately at steady state
  p53 <- p53[p53$time > (max(p53$time) - time_from_end_of_experiment),]
  p21 <- p21[p21$time > (max(p21$time) - time_from_end_of_experiment),]
  
  common_cells <- intersect(unique(p21$variable), unique(p53$variable))
  
  p53_common_cells <- p53[p53$variable %in% common_cells,]
  p21_common_cells <- p21[p21$variable %in% common_cells,]
  # For IR_Nutlin dataset p21 is imaged at a subset of the timepoints that p53 is imaged at
  p53_common_cells <- p53_common_cells[p53_common_cells$time %in% unique(p21_common_cells$time),]
  
  cell <- array(dim = length(unique(p53_common_cells$variable)))
  condition <- array(dim = length(unique(p53_common_cells$variable)))
  p53_mean <- array(dim = length(unique(p53_common_cells$variable)))
  p21_mean <- array(dim = length(unique(p21_common_cells$variable)))
  
  for (cell_number in 1:length(unique(p53_common_cells$variable))) {
    
    cell[cell_number] <- as.character(unique(p21_common_cells$variable)[cell_number])
    
    p53_time <- p53_common_cells[p53_common_cells$variable == cell[cell_number], 'time']
    p53_fluorescence <- p53_common_cells[p53_common_cells$variable == cell[cell_number], 'value']
    p53_mean[cell_number] <- mean(p53_fluorescence)
    
    p21_time <- p21_common_cells[p21_common_cells$variable == cell[cell_number], 'time']
    p21_fluorescence <- p21_common_cells[p21_common_cells$variable == cell[cell_number], 'value']
    p21_mean[cell_number] <- mean(p21_fluorescence)
    
    condition[cell_number] <- p53_common_cells[p53_common_cells$variable == cell[cell_number], 'condition'][1]
  }
  
  integrated_signals <- data.frame(p53_mean, p21_mean, condition, cell)
  
  if (experiment == 'IR_Nutlin') {

    # Fit 2 component Gaussian mixture models to the time averaged p21 values to account for the fact that there are two
    # different p21 steady state levels for each condition due to the two different p21 degradation rates
    mixture_model_ir <- 
      Mclust(log(integrated_signals[integrated_signals$condition == 'IR', 'p21_mean']),
             G = 2, model = 'V', na.rm = TRUE)
    mixture_model_ir_nutlin <- 
      Mclust(log(integrated_signals[integrated_signals$condition == 'IR + nutlin', 'p21_mean']),
             G = 2, model = 'V', na.rm = TRUE)
    
    integrated_signals$mixture <- c(mixture_model_ir$classification, mixture_model_ir_nutlin$classification)
    
    model_ir_1 <- lm(log(p21_mean) ~ log(p53_mean),
                     data = subset(integrated_signals, (condition == 'IR' & mixture == 1)))
    model_ir_2 <- lm(log(p21_mean) ~ log(p53_mean),
                     data = subset(integrated_signals, (condition == 'IR' & mixture == 2)))
    
    model_ir_nutlin_1 <- lm(log(p21_mean) ~ log(p53_mean),
                            data = subset(integrated_signals, (condition == 'IR + nutlin' & mixture == 1)))
    model_ir_nutlin_2 <- lm(log(p21_mean) ~ log(p53_mean),
                            data = subset(integrated_signals, (condition == 'IR + nutlin' & mixture == 2)))
    
    estimates_ir_1 <- summary(model_ir_1)$coefficients[, 1]
    std_errors_ir_1 <- summary(model_ir_1)$coefficients[, 2]
    estimates_ir_2 <- summary(model_ir_2)$coefficients[, 1]
    std_errors_ir_2 <- summary(model_ir_2)$coefficients[, 2]
    
    estimates_ir_nutlin_1 <- summary(model_ir_nutlin_1)$coefficients[, 1]
    std_errors_ir_nutlin_1 <- summary(model_ir_nutlin_1)$coefficients[, 2]
    estimates_ir_nutlin_2 <- summary(model_ir_nutlin_2)$coefficients[, 1]
    std_errors_ir_nutlin_2 <- summary(model_ir_nutlin_2)$coefficients[, 2]
    
    coefficients <- 
      data.frame(
        estimate = c(estimates_ir_1, estimates_ir_2, estimates_ir_nutlin_1, estimates_ir_nutlin_2), 
        std_error = c(std_errors_ir_1, std_errors_ir_2, std_errors_ir_nutlin_1, std_errors_ir_nutlin_2),
        term = rep(c('y-intercept', 'slope'), 4),
        condition = c(rep('IR \n Mixture Component 1', 2), rep('IR \n Mixture Component 2', 2),
                      rep('IR + nutlin \n Mixture Component 1', 2), rep('IR + nutlin \n Mixture Component 2', 2)),
        mixture = c(1, 1, 2, 2, 1, 1, 2, 2),
        colors = c(rep(plot_colors[1], 4), rep(plot_colors[2], 4)))
    
    log_kp_p21_ir_1_estimate <- estimates_ir_1[1] + p21_degradation_rate_mixture_model_fit$mean[2]
    log_kp_p21_ir_2_estimate <- estimates_ir_2[1] + p21_degradation_rate_mixture_model_fit$mean[1]
    log_kp_p21_ir_nutlin_1_estimate <- estimates_ir_nutlin_1[1] + p21_degradation_rate_mixture_model_fit$mean[2]
    log_kp_p21_ir_nutlin_2_estimate <- estimates_ir_nutlin_2[1] + p21_degradation_rate_mixture_model_fit$mean[1]
    
    # Error propagation
    log_kp_p21_ir_1_std_error <- sqrt(std_errors_ir_1[1]^2 + fast_degradation_sigma^2)
    log_kp_p21_ir_2_std_error <- sqrt(std_errors_ir_2[1]^2 + slow_degradation_sigma^2)
    log_kp_p21_ir_nutlin_1_std_error <- sqrt(std_errors_ir_nutlin_1[1]^2 + fast_degradation_sigma^2)
    log_kp_p21_ir_nutlin_2_std_error <- sqrt(std_errors_ir_nutlin_2[1]^2 + slow_degradation_sigma^2)
    
    # Inverse-variance weighted mean
    log_kp_p21_ir_weighted_mean <- 
      (log_kp_p21_ir_1_estimate*(1/log_kp_p21_ir_1_std_error^2) + 
         log_kp_p21_ir_2_estimate*(1/log_kp_p21_ir_2_std_error^2))/
      ((1/log_kp_p21_ir_1_std_error^2) + (1/log_kp_p21_ir_2_std_error^2))
    
    log_kp_p21_ir_nutlin_weighted_mean <- 
      (log_kp_p21_ir_nutlin_1_estimate*(1/log_kp_p21_ir_nutlin_1_std_error^2) + 
         log_kp_p21_ir_nutlin_2_estimate*(1/log_kp_p21_ir_nutlin_2_std_error^2))/
      ((1/log_kp_p21_ir_nutlin_1_std_error^2) + (1/log_kp_p21_ir_nutlin_2_std_error^2))
    
    # Standard error of inverse-variance weighted mean
    log_kp_p21_ir_weighted_sigma <-
      sqrt(1/((1/log_kp_p21_ir_1_std_error^2) + (1/log_kp_p21_ir_2_std_error^2)))
    
    log_kp_p21_ir_nutlin_weighted_sigma <-
      sqrt(1/((1/log_kp_p21_ir_nutlin_1_std_error^2) + (1/log_kp_p21_ir_nutlin_2_std_error^2)))
    
    production_rates <- 
      data.frame(
        estimate = 
          c(log_kp_p21_ir_1_estimate, log_kp_p21_ir_2_estimate, 
            log_kp_p21_ir_nutlin_1_estimate, log_kp_p21_ir_nutlin_2_estimate,
            log_kp_p21_ir_weighted_mean, log_kp_p21_ir_nutlin_weighted_mean),
        std_error = 
          c(log_kp_p21_ir_1_std_error, log_kp_p21_ir_2_std_error,
            log_kp_p21_ir_nutlin_1_std_error, log_kp_p21_ir_nutlin_2_std_error,
            log_kp_p21_ir_weighted_sigma, log_kp_p21_ir_nutlin_weighted_sigma),
        condition =
          c('IR \n Fast degradation', 'IR \n Slow degradation',
            'IR + nutlin \n Fast degradation', 'IR + nutlin \n Slow degradation',
            'IR \n Weighted mean', 'IR + nutlin \n Weighted mean'),
        mixture = c(1, 2, 1, 2, 3, 3),
        colors = c(plot_colors[1], plot_colors[1], plot_colors[2], plot_colors[2], plot_colors[1], plot_colors[2]))
    
  } else if (experiment == 'IR_15min_resolution') {
    
    # Fit 2 component Gaussian mixture models to the time averaged p21 values to account for the fact that there are two
    # different p21 steady state levels for each condition due to the two different p21 degradation rates
    plot_data <- data.frame(
      log_p53 = log(integrated_signals$p53_mean),
      log_p21_to_p53_ratio = log(integrated_signals$p21_mean/integrated_signals$p53_mean))
    
    test_mixture_model <- Mclust(log(integrated_signals[, 'p21_mean']), G = 1:2, model = 'V', na.rm = TRUE)
    print(test_mixture_model$BIC)
    mixture_model <- Mclust(log(integrated_signals[, 'p21_mean']), G = 2, model = 'V', na.rm = TRUE)
    
    integrated_signals$mixture <- mixture_model$classification
    integrated_signals$condition <- 'IR'
    
    regression_model_component_1 <-
      lm(log(p21_mean) ~ log(p53_mean), data = subset(integrated_signals, mixture == 1))
    regression_model_component_2 <-
      lm(log(p21_mean) ~ log(p53_mean), data = subset(integrated_signals, mixture == 2))
    
    estimates_component_1 <- summary(regression_model_component_1)$coefficients[, 1]
    std_errors_component_1 <- summary(regression_model_component_1)$coefficients[, 2]
    estimates_component_2 <- summary(regression_model_component_2)$coefficients[, 1]
    std_errors_component_2 <- summary(regression_model_component_2)$coefficients[, 2]
    
    coefficients <- 
      data.frame(
        estimate = c(estimates_component_1, estimates_component_2), 
        std_error = c(std_errors_component_1, std_errors_component_2),
        term = rep(c('y-intercept', 'slope'), 2),
        condition = c(rep('IR \n Mixture Component 1', 2), rep('IR \n Mixture Component 2', 2)),
        mixture = c(1, 1, 2, 2),
        colors = rep(plot_colors[1], 4))
    
    log_kp_p21_1 <- estimates_component_1[1] + p21_degradation_rate_mixture_model_fit$mean[2]
    log_kp_p21_2 <- estimates_component_2[1] + p21_degradation_rate_mixture_model_fit$mean[1]
    
    production_rates <- data.frame(log_kp_p21_1, log_kp_p21_2)
  }
  
  p_p53_p21 <- 
    ggplot(integrated_signals,
           aes(log(p53_mean), log(p21_mean), 
               color = factor(condition),
               fill = factor(condition),
               shape = factor(mixture))) +
    geom_point() +
    stat_density2d(alpha = 0.5) + 
    geom_smooth(method = 'lm', formula = y ~ x) +
    xlab(TeX('log $\\left( \\bar{ p53(t) } \\right)$')) +
    ylab(TeX('log $\\left( \\bar{ p21(t) } \\right)$')) +
    custom_plot_theme +
    theme(legend.position = c(0.15, 0.8))
  
  p_p21_fits <- 
    ggplot(coefficients,
           aes(x = estimate,
               y = factor(condition),
               color = colors,
               shape = factor(mixture))) +
    geom_point(size = 4, position = position_dodgev(height = 0.5)) + 
    geom_errorbarh(aes(xmin = estimate - std_error, xmax = estimate + std_error),
                   height = 0.2, size = 1, position = position_dodgev(height = 0.5)) +
    facet_grid(~ factor(term, levels = c('slope', 'y-intercept')), scales = 'free') +
    #xlim(c(-1.5, 1.5)) +
    geom_vline(data = subset(coefficients, term == 'slope'), aes(xintercept = 0), linetype = 'dotted') + 
    xlab('Estimate') + ylab('') +
    custom_plot_theme +
    theme(legend.position = 'none')
  
  p_production_rates <- 
    ggplot(production_rates,
           aes(x = estimate,
               y = factor(condition),
               color = colors,
               shape = factor(mixture))) +
    geom_point(
      size = c(4, 4, 8, 4, 4, 8), position = position_dodgev(height = 0.5), alpha = c(0.5, 0.5, 1, 0.5, 0.5, 1)) + 
    geom_errorbarh(
      aes(xmin = estimate - std_error, xmax = estimate + std_error), height = 0.2, size = 1,
      position = position_dodgev(height = 0.5), alpha = c(0.5, 0.5, 1, 0.5, 0.5, 1)) +
    xlab('log(Production rate) (a.u./h)') + ylab('') +
    custom_plot_theme +
    theme(legend.position = 'none')
  
  iffl_abrogation_prediction <- integrated_signals[integrated_signals$condition == 'IR' & integrated_signals$mixture == 1,]
  iffl_abrogation_prediction$condition <- 'Predicted'
  iffl_abrogation_prediction$p21_mean <- iffl_abrogation_prediction$p21_mean*9.386496
  iffl_abrogation_prediction_plot_data <- rbind(integrated_signals, iffl_abrogation_prediction)
  
  p_iffl_abrogation_prediction <- 
    ggplot(iffl_abrogation_prediction_plot_data[
      iffl_abrogation_prediction_plot_data$condition == 'IR' | 
        iffl_abrogation_prediction_plot_data$condition == 'Predicted',]) +
    stat_density2d(data = iffl_abrogation_prediction_plot_data[
      iffl_abrogation_prediction_plot_data$condition == 'IR',],
                   aes(log(p53_mean), log(p21_mean), shape = factor(mixture)),
                   alpha = 0.5, color = plot_colors[1]) +
    geom_point(data = iffl_abrogation_prediction_plot_data[
      iffl_abrogation_prediction_plot_data$condition == 'Predicted',],
               aes(log(p53_mean), log(p21_mean)),
               color = plot_colors[2]) +
    geom_errorbar(data = iffl_abrogation_prediction_plot_data[
      iffl_abrogation_prediction_plot_data$condition == 'Predicted',],
                  aes(x = log(p53_mean),
                      ymin = log(p21_mean - 0.05168517*p21_mean),
                      ymax = log(p21_mean + 0.05168517*p21_mean)),
                  color = plot_colors[2]) +
    xlab(TeX('log $\\left( \\bar{ p53(t) } \\right)$')) +
    ylab(TeX('log $\\left( \\bar{ p21(t) } \\right)$')) +
    xlim(c(5, 9)) + ylim(c(5, 10)) +
    custom_plot_theme +
    theme(legend.position = c(0.15, 0.8))
  
  p53_mean_ir <- mean(integrated_signals[integrated_signals$condition == 'IR', 'p53_mean'])
  p53_se_ir <- sd(integrated_signals[integrated_signals$condition == 'IR', 'p53_mean']/
                    sqrt(length(integrated_signals[integrated_signals$condition == 'IR', 'p53_mean'])))
  p53_mean_ir_nutlin <- mean(integrated_signals[integrated_signals$condition == 'IR + nutlin', 'p53_mean'])
  p53_se_ir_nutlin <- sd(integrated_signals[integrated_signals$condition == 'IR + nutlin', 'p53_mean']/
                           sqrt(length(integrated_signals[integrated_signals$condition == 'IR + nutlin', 'p53_mean'])))
  p53_mean_difference <- p53_mean_ir_nutlin - p53_mean_ir
  p53_se_difference <- sqrt(p53_se_ir^2 + p53_se_ir_nutlin^2)
  
  degradation_abrogation_prediction <- 
    integrated_signals[integrated_signals$condition == 'IR' & integrated_signals$mixture == 1,]
  degradation_abrogation_prediction$condition <- 'Predicted'
  degradation_abrogation_prediction$p53_mean <- degradation_abrogation_prediction$p53_mean + p53_mean_difference
  degradation_abrogation_prediction_plot_data <- rbind(integrated_signals, degradation_abrogation_prediction)
  
  p_degradation_abrogation_prediction <- 
    ggplot(degradation_abrogation_prediction_plot_data[
      degradation_abrogation_prediction_plot_data$condition == 'IR' | 
        degradation_abrogation_prediction_plot_data$condition == 'Predicted',]) +
    stat_density2d(data = degradation_abrogation_prediction_plot_data[
      degradation_abrogation_prediction_plot_data$condition == 'IR',],
      aes(log(p53_mean), log(p21_mean), shape = factor(mixture)),
      alpha = 0.5, color = plot_colors[1]) +
    geom_point(data = degradation_abrogation_prediction_plot_data[
      degradation_abrogation_prediction_plot_data$condition == 'Predicted',],
      aes(log(p53_mean), log(p21_mean)),
      color = plot_colors[2]) +
    geom_errorbarh(data = degradation_abrogation_prediction_plot_data[
      degradation_abrogation_prediction_plot_data$condition == 'Predicted',],
      aes(xmin = log(p53_mean - p53_se_difference),
          xmax = log(p53_mean + p53_se_difference),
          y = log(p21_mean)),
      color = plot_colors[2]) +
    xlab(TeX('log $\\left( \\bar{ p53(t) } \\right)$')) +
    ylab(TeX('log $\\left( \\bar{ p21(t) } \\right)$')) +
    xlim(c(5, 11)) + ylim(c(5, 10)) +
    custom_plot_theme +
    theme(legend.position = c(0.15, 0.8))
  
  print(production_rates)
  
  return(list(p_p53_p21, p_p21_fits, p_production_rates, p_iffl_abrogation_prediction, p_degradation_abrogation_prediction))
}
