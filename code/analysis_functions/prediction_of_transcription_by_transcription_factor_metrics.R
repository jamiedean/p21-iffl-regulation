library(cobs)
library(depmixS4)
library(ggpubr)
library(rms)

source('code/analysis_functions/custom_plot_theme.R')
source('code/analysis_functions/preprocessing.R')

################################################################################################################

smooth_transcription_fractor_dynamics <- function(tf_trace, degree = 2, span = 0.1, plot = FALSE) {
  
  if (plot == TRUE) {
    
    tf_trace$smooth <- loess(value ~ time, degree = degree, span = span, data = tf_trace)$fitted
    
    p <- ggplot(tf_trace, aes(time, value)) +
      geom_line(linewidth = 1.5, alpha = 0.5, color = 'gray') +
      geom_line(aes(time, smooth), linewidth = 1.5, color = plot_colors[1]) +
      xlab('Time (h)') + ylab('p53 (a.u.)') +
      custom_plot_theme
    
    return(p)
    
  } else {
    
    tf_trace$value <- loess(value ~ time, degree = degree, span = span, data = tf_trace)$fitted
   
    return(tf_trace) 
  }
}

detrend_transcription_factor_dynamics <- function(tf_trace, plot = FALSE, num_metrics = 2) {
  
  bspline_number_knots <- 4
  bspline_degree <- 2
  
  bspline_regression <- 
    cobs(tf_trace$time, tf_trace$value, nknots = bspline_number_knots, degree = bspline_degree, knots.add = TRUE)
  tf_trace$trend <- bspline_regression$fitted
  tf_trace$detrended <- tf_trace$value - tf_trace$trend
  
  if (plot == TRUE) {
    
    p_trend <- ggplot(tf_trace, aes(time, trend)) +
      geom_line(aes(time, value), linewidth = 1.5, color = plot_colors[1], alpha = 0.5) +
      geom_line(linewidth = 1.5, color = 'black') +
      xlab('Time (h)') + ylab('p53 (a.u.)') +
      custom_plot_theme
    
    if (num_metrics == 2) {
      
      line_color <- plot_colors[1]
      
    } else if (num_metrics == 3) {
      
      line_color <- plot_colors[2]
    }
    
    p_detrended <- ggplot(tf_trace, aes(time, detrended)) +
      geom_line(linewidth = 1.5, color = line_color) +
      xlab('Time (h)') + ylab('p53 detrended (a.u.)') +
      custom_plot_theme
    
    return(list(p_trend, p_detrended))
    
  } else {
    
    tf_detrended_df <- data.frame(time = tf_trace$time, value = tf_trace$detrended)
    
    return(tf_detrended_df) 
  }
}

infer_gene_state <- function(ms2_trace, plot = FALSE, num_metrics = 2) {
  
  hidden_markov_model <- depmix(value ~ 1, data = ms2_trace, nstates = 2, family = gaussian())
  fit_hidden_markov_model <- fit(hidden_markov_model)
  
  # Predict the states by estimating the posterior
  ms2_trace$estimated_states <- posterior(fit_hidden_markov_model, type = 'viterbi')
  
  if (mean(ms2_trace$value[which(ms2_trace$estimated_states$state == 1)]) > 
      mean(ms2_trace$value[which(ms2_trace$estimated_states$state == 2)])) {
    
    ms2_trace$estimated_states$state[ms2_trace$estimated_states$state == 1] <- 1
    ms2_trace$estimated_states$state[ms2_trace$estimated_states$state == 2] <- 0
    
  } else {
    
    ms2_trace$estimated_states$state[ms2_trace$estimated_states$state == 1] <- 0
    ms2_trace$estimated_states$state[ms2_trace$estimated_states$state == 2] <- 1
  }
  
  if (plot == TRUE) {
    
    if (num_metrics == 2) {
      
      line_color <- plot_colors[3]
      
    } else if (num_metrics == 3) {
      
      line_color <- plot_colors[5]
    }
    
    p <- ggplot(ms2_trace, aes(time, value)) +
      geom_line(linewidth = 1.5, alpha = 0.5, color = 'gray') +
      geom_line(aes(time, ms2_trace$estimated_states$state), linewidth = 1.5, color = line_color) +
      xlab('Time (h)') + ylab('p21-MS2 (a.u.)') +
      custom_plot_theme
    
    return(p)
    
  } else {
    
    return(ms2_trace) 
  }
}

calculate_change_in_tf_level <- function(tf_cell, plot = FALSE, num_metrics = 2, increase_only = TRUE) {
  
  tf_cell_change <- diff(tf_cell$value)
  
  if (increase_only == TRUE) {
    
    tf_cell_change[tf_cell_change < 0] <- 0
    ylab <- 'p53 increase (a.u.)'
    
    if (num_metrics == 2) {
      
      line_color <- plot_colors[3]
      
    } else if (num_metrics == 3) {
      
      line_color <- plot_colors[4]
    }
    
  } else if (increase_only == FALSE) {
    
    ylab <- 'p53 change (a.u.)'
    
    if (num_metrics == 2) {
      
      line_color <- plot_colors[2]
      
    } else if (num_metrics == 3) {
      
      line_color <- plot_colors[3]
    }
  }
  
  tf_cell_change_df <- data.frame(time = tf_cell$time[-1], value = tf_cell_change)
  
  if (plot == TRUE) {
    
    p <- ggplot(tf_cell_change_df, aes(time, tf_cell_change)) +
      geom_line(linewidth = 1.5, color = line_color) +
      xlab('Time (h)') + ylab(ylab) +
      custom_plot_theme
    
    return(p)
    
  } else {
   
    return(tf_cell_change_df) 
  }
}

shift_traces <- function(traces_data_cell) {
  
  # Maximum shift cannot be greater than half the p53 oscillation period (5.5 h in humans)
  maximum_lag <- (5.5/2)/(traces_data_cell$time[2] - traces_data_cell$time[1])
  
  cross_correlation <- ccf(traces_data_cell$tf_metric, traces_data_cell$gene_status, lag.max = maximum_lag)
  optimal_lag <- cross_correlation$lag[which.max(cross_correlation$acf)]
  
  if (optimal_lag <= 0) {
    
    traces_data_cell_shifted <- 
      data.frame(
        tf_metric = traces_data_cell$tf_metric[1:(length(traces_data_cell$tf_metric) + optimal_lag)],
        gene_status = traces_data_cell$gene_status[(1 - optimal_lag):length(traces_data_cell$tf_metric)])
  } else if (optimal_lag > 0) {
    
    traces_data_cell_shifted <- 
      data.frame(
        tf_metric = traces_data_cell$tf_metric[1:(length(traces_data_cell$tf_metric) - optimal_lag)],
        gene_status = 
          traces_data_cell$gene_status[optimal_lag:(length(traces_data_cell$tf_metric) - 1)])
  }
  
  return(traces_data_cell_shifted)
}

fit_logistic_regression_model <- function(traces_data_shifted, tf_metric_name, num_metrics = 2, plot = FALSE) {
  
  data_dist <<- datadist(traces_data_shifted)
  options(datadist = 'data_dist')
  
  # Small positive penalty term added to avoid singular matrix error (needed for nutlin dataset)
  logistic_regression <- lrm(gene_status ~ tf_metric, data = traces_data_shifted, maxit = 100, penalty = 1e-1)
  print(logistic_regression)
  auc <- logistic_regression[[3]]['C']
  coefficient <- coef(logistic_regression)[['tf_metric']]
  
  if (plot == TRUE) {
    
    if (num_metrics == 2) {
      
      if (tf_metric_name == 'p53 level') {
        metric_dependent_color <- plot_colors[1]
      } else if (tf_metric_name == 'p53 change') {
        metric_dependent_color <- plot_colors[2]
      }
    } else if (num_metrics == 3) {
      
      if (tf_metric_name == 'p53 level') {
        metric_dependent_color <- plot_colors[1]
      } else if (tf_metric_name == 'p53 detrended') {
        metric_dependent_color <- plot_colors[2]
      } else if (tf_metric_name == 'p53 change') {
        metric_dependent_color <- plot_colors[3]
      } else if (tf_metric_name == 'p53 increase') {
        metric_dependent_color <- plot_colors[4]
      }
    }
    
    p <- ggplot(traces_data_shifted, aes(x = tf_metric, y = as.numeric(gene_status) - 1)) + 
      geom_rug(data = traces_data_shifted[(as.numeric(traces_data_shifted$gene_status) - 1) == 0,],
               alpha = 0.3, color = metric_dependent_color, sides = 'b') +
      geom_rug(data = traces_data_shifted[(as.numeric(traces_data_shifted$gene_status) - 1) == 1,],
               alpha = 0.3, color = metric_dependent_color, sides = 't') +
      geom_histogram(data = traces_data_shifted[(as.numeric(traces_data_shifted$gene_status) - 1) == 0,],
                     aes(x = tf_metric, y = after_stat(count)/100), bins = 25, alpha = 0.3, color = 'white',
                     fill = metric_dependent_color) +
      geom_histogram(data = traces_data_shifted[(as.numeric(traces_data_shifted$gene_status) - 1) == 1,],
                     aes(x = tf_metric, y = -after_stat(count)/100), bins = 25, alpha = 0.3, color = 'white',
                     fill = metric_dependent_color, position = position_nudge(y = 1)) +
      stat_smooth(method = 'glm', se = FALSE, method.args = list(family = binomial), linewidth = 1.5,
                  color = metric_dependent_color) + 
      annotate('text', label = paste('AUC = ', round(auc, 2)), size = 5,
               x = 0.7*max(traces_data_shifted$tf_metric), y = 0.25) +
      xlab(paste(tf_metric_name, '(a.u.)')) + ylab('p21 gene status') +
      custom_plot_theme

    return(p)
    
  } else {
      
    return(auc)
  }
}

predict_promoter_status <- function(tf, ms2, cell_number, tf_metric = 'change', num_metrics = 2, detrend_first = TRUE,
                                    plot = FALSE) {
  
  cells <- intersect(unique(tf$variable), unique(ms2$variable))
  
  tf_cell <- tf[tf$variable == cells[cell_number], c('time', 'value')]
  ms2_cell <- ms2[ms2$variable == cells[cell_number], c('time', 'value')]
  
  tf_cell <- smooth_transcription_fractor_dynamics(tf_cell, plot = FALSE)
  tf_cell_detrended <- detrend_transcription_factor_dynamics(tf_cell, plot = FALSE)
  
  if (detrend_first == TRUE) {
    
    tf_cell_change <- calculate_change_in_tf_level(tf_cell_detrended, num_metrics = num_metrics, increase_only = FALSE)
    tf_cell_increase <- calculate_change_in_tf_level(tf_cell_detrended, num_metrics = num_metrics, increase_only = TRUE)
  
  } else if (detrend_first == FALSE) {
    
    tf_cell_change <- calculate_change_in_tf_level(tf_cell, num_metrics = num_metrics, increase_only = FALSE)
    tf_cell_increase <- calculate_change_in_tf_level(tf_cell, num_metrics = num_metrics, increase_only = TRUE)
  }

  ms2_cell <- infer_gene_state(ms2_cell, plot = FALSE, num_metrics = num_metrics)

  if (tf_metric == 'increase') {
    
    traces_data <- 
      data.frame(time = tf_cell$time[-1],
                 tf_metric = tf_cell_increase$value,
                 gene_status = factor(ms2_cell$estimated_states$state[-1]))
    tf_metric_name <- 'p53 increase'
    
  } else if (tf_metric == 'level') {
    
    if (detrend_first == TRUE) {
      
      traces_data <- 
        data.frame(time = tf_cell$time,
                   tf_metric = tf_cell_detrended$value,
                   gene_status = factor(ms2_cell$estimated_states$state))
      
    } else if (detrend_first == FALSE) {
      
      traces_data <- 
        data.frame(time = tf_cell$time,
                   tf_metric = tf_cell$value,
                   gene_status = factor(ms2_cell$estimated_states$state))
    }
    
    tf_metric_name <- 'p53 level'
    
  } else if (tf_metric == 'change') {
    
    traces_data <- 
      data.frame(time = tf_cell$time[-1],
                 tf_metric = tf_cell_change$value,
                 gene_status = factor(ms2_cell$estimated_states$state[-1]))
    tf_metric_name <- 'p53 change'
    
  } else if (tf_metric == 'detrended') {
    
    traces_data <- 
      data.frame(time = tf_cell$time,
                 tf_metric = tf_cell_detrended$value,
                 gene_status = factor(ms2_cell$estimated_states$state))
    tf_metric_name <- 'p53 detrended'
  }

  traces_data_shifted <- shift_traces(traces_data)
  
  if (length(unique(traces_data_shifted$gene_status)) == 1) {
    
    auc <- NA
    
  } else {
    
    auc <- fit_logistic_regression_model(traces_data_shifted, tf_metric_name, num_metrics = num_metrics, plot = plot)
  }
  
  return(auc)
}

compare_predictive_performance_of_metrics_of_transcription <- function(tf, ms2, num_metrics = 2, detrend_first = TRUE) {
  
  cells <- intersect(unique(tf$variable), unique(ms2$variable))
  auc_level <- array(dim = length(cells))
  auc_detrended <- array(dim = length(cells))
  auc_change <- array(dim = length(cells))
  
  for (cell_number in 1:length(cells)) {
    
    print(cell_number)
    if (detrend_first == FALSE) {
      
      auc_level[cell_number] <- 
        predict_promoter_status(
          tf, ms2, cell_number, tf_metric = 'level', num_metrics = num_metrics, detrend_first = detrend_first,
          plot = FALSE)
    }
    
    if (num_metrics == 3 | (num_metrics == 2 & detrend_first == TRUE)) {
      
      auc_detrended[cell_number] <- 
        predict_promoter_status(
          tf, ms2, cell_number, tf_metric = 'detrended', num_metrics = num_metrics, detrend_first = detrend_first,
          plot = FALSE)
    }
    
    auc_change[cell_number] <- 
      predict_promoter_status(
        tf, ms2, cell_number, tf_metric = 'change', num_metrics = num_metrics, detrend_first = detrend_first,
        plot = FALSE)
  }
    
  if (num_metrics == 2) {
    
    if (detrend_first == TRUE) {
      
      model_comparison <- 
        data.frame(
          cell = cells,
          metric = factor(c(rep('p53 level', length(auc_detrended)),
                            rep('p53 change', length(auc_change))),
                          levels = c('p53 level', 'p53 change')),
          auc = c(auc_detrended, auc_change))
      
      print(mean(auc_detrended, na.rm = TRUE))
      print(sd(auc_detrended, na.rm = TRUE))
      
    } else if (detrend_first == FALSE) {
      
      model_comparison <- 
        data.frame(
          cell = cells,
          metric = factor(c(rep('p53 level', length(auc_level)),
                            rep('p53 change', length(auc_change))),
                          levels = c('p53 level', 'p53 change')),
          auc = c(auc_level, auc_change))
      
      print(mean(auc_level, na.rm = TRUE))
      print(sd(auc_level, na.rm = TRUE))
    }
    
    print(mean(auc_change))
    print(sd(auc_change))
    
  } else if (num_metrics == 3) {
    
    print(mean(auc_level, na.rm = TRUE))
    print(sd(auc_level, na.rm = TRUE))
    print(mean(auc_detrended, na.rm = TRUE))
    print(sd(auc_detrended, na.rm = TRUE))
    print(mean(auc_change))
    print(sd(auc_change))
    
    model_comparison <- 
      data.frame(
        cell = cells,
        metric = factor(c(rep('p53 level', length(auc_level)),
                          rep('p53 detrended', length(auc_detrended)),
                          rep('p53 change', length(auc_change))),
                        levels = c('p53 level', 'p53 detrended', 'p53 change')),
        auc = c(auc_level, auc_detrended, auc_change))
  }
  
  cells_with_na <- model_comparison[which(is.na(model_comparison$auc)), 'cell']

  if (length(cells_with_na > 0)) {
    model_comparison <- model_comparison[model_comparison$cell != cells_with_na,] 
  }
  
  p <- ggplot(model_comparison, aes(x = metric, y = auc)) +
    geom_violin(aes(fill = metric)) +
    geom_boxplot(width = 0.35) +
    xlab('Metric') + ylab('AUC') +
    custom_plot_theme + theme(legend.position = 'none')
  
  if (num_metrics == 2) {
      
    my_comparisons <- NULL
    label_x = 1.25
    label_y <- c(1.1)
    y_axis_limits <- c(0.25, 1.2)
    
  } else if (num_metrics == 3) {
      
    my_comparisons <- list(c('p53 level', 'p53 detrended'),
                           c('p53 detrended', 'p53 change'),
                           c('p53 level', 'p53 change')) 
    label_x = 1.25
    label_y <- c(1.0, 1.1, 1.2)
    y_axis_limits <- c(0.25, 1.3)
  }
    
  p <- p + 
    ylim(y_axis_limits) +
    stat_compare_means(
      comparison = my_comparisons, paired = TRUE, method = 'wilcox.test', label.x = label_x, label.y = label_y, size = 5)
  
  return(p)
}
