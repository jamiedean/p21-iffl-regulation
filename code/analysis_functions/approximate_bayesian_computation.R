library(abcrf)
library(parallel)
library(scales)

source('code/analysis_functions/autocorrelation_analysis.R')
source('code/analysis_functions/preprocessing.R')
source('code/stochastic_models/p21rna_ssa.R')

################################################################################################################

generate_correlation_functions <- function(p53_pred, p21rna_pred, p21_pred) {
  
  p53_autocorrelation_pred <- NULL
  p21rna_autocorrelation_pred <- NULL
  p21_autocorrelation_pred <- NULL
  p53_p21rna_crosscorrelation_pred <- NULL
  p21rna_p21_crosscorrelation_pred <- NULL
  p53_p21_crosscorrelation_pred <- NULL
  
  if (is.null(p53_pred) == FALSE) {
    p53_autocorrelation_pred <- compute_correlation_function(p53_pred, p53_pred)
  }
  
  if (is.null(p21rna_pred) == FALSE) {
    p21rna_autocorrelation_pred <- compute_correlation_function(p21rna_pred, p21rna_pred)
  }
  
  if (is.null(p21_pred) == FALSE) {
    p21_autocorrelation_pred <- compute_correlation_function(p21_pred, p21_pred)
  }
  
  if (is.null(p53_pred) == FALSE & is.null(p21rna_pred) == FALSE) {
    p53_p21rna_crosscorrelation_pred <- compute_correlation_function(p53_pred, p21rna_pred)
  }
  
  if (is.null(p21rna_pred) == FALSE & is.null(p21_pred) == FALSE) {
    p21rna_p21_crosscorrelation_pred <- compute_correlation_function(p21rna_pred, p21_pred)
  }
  
  if (is.null(p53_pred) == FALSE & is.null(p21_pred) == FALSE) {
    p53_p21_crosscorrelation_pred <- compute_correlation_function(p53_pred, p21_pred)
  }
  
  return(list(p53_autocorrelation_pred, p21rna_autocorrelation_pred, p21_autocorrelation_pred,
              p53_p21rna_crosscorrelation_pred, p21rna_p21_crosscorrelation_pred, p53_p21_crosscorrelation_pred))
}

generate_reference_table <- function(experiment_pred, model, parameter_sets, model_type, plot = FALSE) {
  
  num_cores <- detectCores()
  
  # If running on laptop do not use all available cores
  if (num_cores <= 8) {
    cluster <- makeCluster(num_cores - 2, type = 'FORK')
  } else if (num_cores > 8) {
    cluster <- makeCluster(num_cores, type = 'FORK')
  }
  
  simulated_datasets <- parSapply(cluster, 1:dim(parameter_sets)[1], function(parameter_set) {
    generate_simulated_datasets(experiment_pred, model, parameter_sets[parameter_set,], model_type)
  })
  
  clusterExport(cluster, 'simulated_datasets', envir = environment())
  
  mean_sd_timecourses <- parSapply(cluster, 1:dim(parameter_sets)[1], function(parameter_set) {
    generate_mean_sd_timecourses(simulated_datasets[1, parameter_set][[1]],
                                 simulated_datasets[2, parameter_set][[1]],
                                 simulated_datasets[3, parameter_set][[1]])
  })
  
  correlation_functions <- parSapply(cluster, 1:dim(parameter_sets)[1], function(parameter_set) {
    generate_correlation_functions(simulated_datasets[1, parameter_set][[1]],
                                   simulated_datasets[2, parameter_set][[1]],
                                   simulated_datasets[3, parameter_set][[1]])
  })
  
  stopCluster(cluster)
  
  p53_means <- matrix(nrow = dim(parameter_sets)[1], ncol = length(unique(simulated_datasets[[1]]$time)))
  p21rna_means <- matrix(nrow = dim(parameter_sets)[1], ncol = length(unique(simulated_datasets[[2]]$time)))
  p21_means <- matrix(nrow = dim(parameter_sets)[1], ncol = length(unique(simulated_datasets[[3]]$time)))
  p53_sds <- matrix(nrow = dim(parameter_sets)[1], ncol = length(unique(simulated_datasets[[1]]$time)))
  p21rna_sds <- matrix(nrow = dim(parameter_sets)[1], ncol = length(unique(simulated_datasets[[2]]$time)))
  p21_sds <- matrix(nrow = dim(parameter_sets)[1], ncol = length(unique(simulated_datasets[[3]]$time)))
  p53_correlations <- 
    matrix(nrow = dim(parameter_sets)[1], ncol = length(unique(simulated_datasets[[1]]$time)))
  p21rna_correlations <- 
    matrix(nrow = dim(parameter_sets)[1], ncol = length(unique(simulated_datasets[[2]]$time)))
  p21_correlations <- 
    matrix(nrow = dim(parameter_sets)[1], ncol = length(unique(simulated_datasets[[3]]$time)))
  
  for (parameter_set in 1:dim(parameter_sets)[1]) {
    
    p53_means[parameter_set,] <- mean_sd_timecourses[1, parameter_set][[1]]$mean
    p21rna_means[parameter_set,] <- mean_sd_timecourses[2, parameter_set][[1]]$mean
    p21_means[parameter_set,] <- mean_sd_timecourses[3, parameter_set][[1]]$mean
    p53_sds[parameter_set,] <- mean_sd_timecourses[1, parameter_set][[1]]$sd
    p21rna_sds[parameter_set,] <- mean_sd_timecourses[2, parameter_set][[1]]$sd
    p21_sds[parameter_set,] <- mean_sd_timecourses[3, parameter_set][[1]]$sd
    p53_correlations[parameter_set,] <- correlation_functions[1, parameter_set][[1]]$G_tau
    p21rna_correlations[parameter_set,] <- correlation_functions[2, parameter_set][[1]]$G_tau
    p21_correlations[parameter_set,] <- correlation_functions[3, parameter_set][[1]]$G_tau
  }
  
  write.table(parameter_sets,
              paste(paste(paste('outputs/', experiment_pred, sep = ''),
                          model, sep = '_'), '_parameter_sets.csv', sep = ''),
              sep = ',', col.names = FALSE, row.names = FALSE)
  
  write.table(p53_means,
              paste(paste(paste('outputs/', experiment_pred, sep = ''),
                          model, sep = '_'), '_p53_mean_timecourses.csv', sep = ''),
              sep = ',', col.names = FALSE, row.names = FALSE)
  
  write.table(p21rna_means,
              paste(paste(paste('outputs/', experiment_pred, sep = ''),
                          model, sep = '_'), '_p21rna_mean_timecourses.csv', sep = ''),
              sep = ',', col.names = FALSE, row.names = FALSE)
  
  write.table(p21_means,
              paste(paste(paste('outputs/', experiment_pred, sep = ''),
                          model, sep = '_'), '_p21_mean_timecourses.csv', sep = ''),
              sep = ',', col.names = FALSE, row.names = FALSE)
  
  write.table(p53_sds,
              paste(paste(paste('outputs/', experiment_pred, sep = ''),
                          model, sep = '_'), '_p53_sd_timecourses.csv', sep = ''),
              sep = ',', col.names = FALSE, row.names = FALSE)
  
  write.table(p21rna_sds,
              paste(paste(paste('outputs/', experiment_pred, sep = ''),
                          model, sep = '_'), '_p21rna_sd_timecourses.csv', sep = ''),
              sep = ',', col.names = FALSE, row.names = FALSE)
  
  write.table(p21_sds,
              paste(paste(paste('outputs/', experiment_pred, sep = ''),
                          model, sep = '_'), '_p21_sd_timecourses.csv', sep = ''),
              sep = ',', col.names = FALSE, row.names = FALSE)
  
  write.table(p53_correlations,
              paste(paste(paste('outputs/', experiment_pred, sep = ''),
                          model, sep = '_'), '_p53_correlation_functions.csv', sep = ''),
              sep = ',', col.names = FALSE, row.names = FALSE)
  
  write.table(p21rna_correlations,
              paste(paste(paste('outputs/', experiment_pred, sep = ''),
                          model, sep = '_'), '_p21rna_correlation_functions.csv', sep = ''),
              sep = ',', col.names = FALSE, row.names = FALSE)
  
  write.table(p21_correlations,
              paste(paste(paste('outputs/', experiment_pred, sep = ''),
                          model, sep = '_'), '_p21_correlation_functions.csv', sep = ''),
              sep = ',', col.names = FALSE, row.names = FALSE)
}

compute_summary_statistics <- function(p53_autocorrelation = NULL, p21rna_autocorrelation = NULL,
                                       p21_autocorrelation = NULL, p53_p21rna_crosscorrelation = NULL,
                                       p21rna_p21_crosscorrelation = NULL, p53_p21_crosscorrelation = NULL) {
  # G intercept, tau intercept
  summary_statistics <- array(dim = 6)
  
  if (is.null(p53_autocorrelation) == FALSE) {
    summary_statistics[1] <- p53_autocorrelation$G_tau[1]
  }
  if (is.null(p21rna_autocorrelation) == FALSE) {
    summary_statistics[2] <- p21rna_autocorrelation$G_tau[1]
  }
  if (is.null(p21_autocorrelation) == FALSE) {
    summary_statistics[3] <- p21_autocorrelation$G_tau[1]
  }
  if (is.null(p53_p21rna_crosscorrelation) == FALSE) {
    summary_statistics[4] <- p53_p21rna_crosscorrelation$G_tau[1]
  }
  if (is.null(p21rna_p21_crosscorrelation) == FALSE) {
    summary_statistics[5] <- p21rna_p21_crosscorrelation$G_tau[1]
  }
  if (is.null(p53_p21_crosscorrelation) == FALSE) {
    summary_statistics[6] <- p53_p21_crosscorrelation$G_tau[1]
  }
  
  return(summary_statistics)
}

save_reference_table <- function(experiment_pred, condition, model1, model2) {
  
  if (experiment_pred == 'Simulated_p53_2min' | 
      experiment_pred == 'Simulated_p53_15min') {
    
    model1_reference_table <- 
      read.table(paste(paste(paste(
        'outputs/', experiment_pred, sep = ''),
        model1, sep = '_'),
        '_p53_correlation_functions.csv', sep = ''), sep = ',')
    model2_reference_table <- 
      read.table(paste(paste(paste(
        'outputs/', experiment_pred, sep = ''),
        model2, sep = '_'),
        '_p53_correlation_functions.csv', sep = ''), sep = ',')
    
  } else if (experiment_pred == 'simulated_p21rna_2min' | 
             experiment_pred == 'IR_15min_resolution' |
             experiment_pred == 'IR_Nutlin') {
    
    n_simulations <- 5000
    
    model1_p21rna_acf <- c()
    model2_p21rna_acf <- c()
    model1_p53_p21rna_ccf <- c()
    model2_p53_p21rna_ccf <- c()
    
    for (file_number in 1:n_simulations) {
      
      model1_p21rna_acf_file <- read.csv(paste(paste(paste(paste(paste(paste(paste(
        'outputs/correlation_functions/', experiment_pred, sep = ''),
        condition, sep = '/'),
        model1, sep = '/'), sep = ''),
        'p21rna_autocorrelation', sep = '/'),
        file_number, sep = '_'), '.csv', sep = ''))
      model1_p21rna_acf <- rbind(model1_p21rna_acf, t(model1_p21rna_acf_file$x))
      
      model1_p53_p21rna_ccf_file <- read.csv(paste(paste(paste(paste(paste(paste(paste(
        'outputs/correlation_functions/', experiment_pred, sep = ''),
        condition, sep = '/'),
        model1, sep = '/'), sep = ''),
        'p53_p21rna_crosscorrelation', sep = '/'),
        file_number, sep = '_'), '.csv', sep = ''))
      model1_p53_p21rna_ccf <- rbind(model1_p53_p21rna_ccf, t(model1_p53_p21rna_ccf_file$x))
      
      model2_p21rna_acf_file <- read.csv(paste(paste(paste(paste(paste(paste(paste(
        'outputs/correlation_functions/', experiment_pred, sep = ''),
        condition, sep = '/'),
        model2, sep = '/'), sep = ''),
        'p21rna_autocorrelation', sep = '/'),
        file_number, sep = '_'), '.csv', sep = ''))
      model2_p21rna_acf <- rbind(model2_p21rna_acf, t(model2_p21rna_acf_file$x))
      
      model2_p53_p21rna_ccf_file <- read.csv(paste(paste(paste(paste(paste(paste(paste(
        'outputs/correlation_functions/', experiment_pred, sep = ''),
        condition, sep = '/'),
        model2, sep = '/'), sep = ''),
        'p53_p21rna_crosscorrelation', sep = '/'),
        file_number, sep = '_'), '.csv', sep = ''))
      model2_p53_p21rna_ccf <- rbind(model2_p53_p21rna_ccf, t(model2_p53_p21rna_ccf_file$x))
    }
    
    colnames(model1_p21rna_acf) <- paste('p21rna_acf_', col(model1_p21rna_acf)[1,], sep = '')
    colnames(model1_p53_p21rna_ccf) <- paste('p53_p21rna_ccf_', col(model1_p53_p21rna_ccf)[1,], sep = '')
    colnames(model2_p21rna_acf) <- paste('p21rna_acf_', col(model2_p21rna_acf)[1,], sep = '')
    colnames(model2_p53_p21rna_ccf) <- paste('p53_p21rna_ccf_', col(model2_p53_p21rna_ccf)[1,], sep = '')

    model1_reference_table <- cbind(model1_p21rna_acf, model1_p53_p21rna_ccf)
    model2_reference_table <- cbind(model2_p21rna_acf, model2_p53_p21rna_ccf)
    combined_models <- data.frame(rbind(model1_reference_table, model2_reference_table))
    combined_models$model_index <-
      as.factor(c(rep(1, dim(model1_reference_table)[1]), rep(2, dim(model2_reference_table)[1])))

    write.csv(combined_models, paste(paste(paste(paste(paste(
      'outputs/reference_tables/', experiment_pred, sep = ''),
      condition, sep = '_'),
      model1, sep = '_'),
      model2, sep = '_'),
      '.csv', sep = ''), row.names = FALSE)
    
  } else if (experiment_pred == 'Simulated_p21rna_p21_2min' | 
             experiment_pred == 'Simulated_p21rna_p21_15min') {
    
    model1_reference_table <- 
      read.table(paste(paste(paste(
        'outputs/', experiment_pred, sep = ''),
        model1, sep = '_'),
        '_p21rna_correlation_functions.csv', sep = ''), sep = ',')
    model2_reference_table <- 
      read.table(paste(paste(paste(
        'outputs/', experiment_pred, sep = ''),
        model2, sep = '_'),
        '_p21rna_correlation_functions.csv', sep = ''), sep = ',')
  }
}

perform_model_selection <- function(experiment_meas, experiment_pred, condition, model1, model2) {
    
  if (experiment_meas == 'IR_Nutlin') {
    
    p53_meas <- load_data(experiment_meas, 'traces_CFP', FALSE)
    p21rna_meas <- load_data(experiment_meas, 'traces_YFP', FALSE)
    p53_meas <- p53_meas[p53_meas$condition == condition,]
    p21rna_meas <- p21rna_meas[p21rna_meas$condition == condition,]
    
  } else {
    
    p53_meas <- load_data(experiment_meas, 'traces_CFP', FALSE)
    p21rna_meas <- load_data(experiment_meas, 'traces_MS2', FALSE)
  }
    
  # Fluroescence maturation correction
  # CFP has ~ 50 min maturation time
  p53_meas <- p53_meas[p53_meas$time > 50/60,]
  p21rna_meas <- p21rna_meas[p21rna_meas$time <= p21rna_meas$time[sum(unique(p21rna_meas$time) >= 50/60)],]
  p21rna_meas$time <- p21rna_meas$time + 50/60
    
  meas_p21rna_acf <- data.frame(t(compute_correlation_function(p21rna_meas, p21rna_meas)))
  meas_p53_p21rna_ccf <- data.frame(t(compute_correlation_function(p53_meas, p21rna_meas)))
  colnames(meas_p21rna_acf) <- gsub('V', 'p21rna_acf_', colnames(meas_p21rna_acf))
  colnames(meas_p53_p21rna_ccf) <- gsub('V', 'p53_p21rna_ccf_', colnames(meas_p53_p21rna_ccf))
  measured_summary_statistics <- cbind(meas_p21rna_acf, meas_p53_p21rna_ccf)
  turning_points <- which(diff(sign(diff(c(t(meas_p21rna_acf))))) != 0)[1:5]
  measured_summary_statistics <- measured_summary_statistics[turning_points]
  
  combined_models <- read.csv(paste(paste(paste(paste(paste(
    'outputs/reference_tables/', experiment_pred, sep = ''),
    condition, sep = '_'),
    model1, sep = '_'),
    model2, sep = '_'),
    '.csv', sep = ''))
  combined_models <- combined_models[, c(turning_points, ncol(combined_models))]
  combined_models$model_index <- as.factor(combined_models$model_index)
  
  model_comparison <- 
    abcrf(model_index ~., data = combined_models, lda = TRUE, ntree = 3000, paral = TRUE)
  print(model_comparison)
  err.abcrf(model_comparison, training = combined_models)

  plot(model_comparison, training = combined_models, obs = measured_summary_statistics)
  model_selection <- predict(model_comparison, obs = measured_summary_statistics, training = combined_models)
  print(model_selection)
  
  return(model_selection)
}

plot_correlation_function_fit <- function(experiment_meas, experiment_pred, condition, model,
                                          correlation_function_name) {
  
  if (experiment_meas == 'IR_Nutlin') {
    
    p53_meas <- load_data(experiment_meas, 'traces_CFP', FALSE)
    p21rna_meas <- load_data(experiment_meas, 'traces_YFP', FALSE)
    
    if (condition %in% c('IR', 'IR + nutlin')) {
      
      p53_meas <- p53_meas[p53_meas$condition == condition,]
      p21rna_meas <- p21rna_meas[p21rna_meas$condition == condition,]
    }
    
  } else if (experiment_meas == 'IR_15min_resolution') {
    
    p53_meas <- load_data(experiment_meas, 'traces_CFP', FALSE)
    p21rna_meas <- load_data(experiment_meas, 'traces_MS2', FALSE)
    
    if (condition %in% c('2.5 Gy', '5 Gy', '10 Gy')) {
     
      p53_meas <- p53_meas[p53_meas$condition == condition,]
      p21rna_meas <- p21rna_meas[p21rna_meas$condition == condition,]
    }
    
  } else {
    
    p53_meas <- load_data(experiment_meas, 'traces_CFP', FALSE)
    p21rna_meas <- load_data(experiment_meas, 'traces_MS2', FALSE) 
  }
  
  # Fluroescence maturation correction
  # CFP has ~ 50 min maturation time
  p53_meas <- p53_meas[p53_meas$time > 50/60,]
  p21rna_meas <- p21rna_meas[p21rna_meas$time <= p21rna_meas$time[sum(unique(p21rna_meas$time) > 50/60)],]
  p21rna_meas$time <- p21rna_meas$time + 50/60
  
  if (model == 'telegraph_standard') {
    
    model_index <- 1
    model_label <- 'PR Model'
    line_color <- plot_colors[1]
    ribbon_color <- plot_colors[1]
    
  } else if (model == 'telegraph_iffl') {
    
    model_index <- 2
    model_label <- 'IFFL Model'
    line_color <- plot_colors[2]
    ribbon_color <- plot_colors[2]
  }
  
  model1 <- 'telegraph_standard'
  model2 <- 'telegraph_iffl'
  
  combined_models <- read.csv(paste(paste(paste(paste(paste(
    'outputs/reference_tables/', experiment_pred, sep = ''),
    condition, sep = '_'),
    model1, sep = '_'),
    model2, sep = '_'),
    '.csv', sep = ''))
  
  reference_table <- 
    combined_models[combined_models$model_index == model_index,
                    !(names(combined_models) %in% c('model_index'))]
  
  sim_correlation_function <- reference_table[, grep(correlation_function_name, names(reference_table))]
  
  if (correlation_function_name == 'p21rna_acf') {
    
    delay_times <- seq(0, (dim(sim_correlation_function)[2] - 1))*(p21rna_meas$time[2] - p21rna_meas$time[1])
    meas_correlation_function <- data.frame(t(compute_correlation_function(p21rna_meas, p21rna_meas)))
    colnames(meas_correlation_function) <- gsub('V', 'p21rna_acf_', colnames(meas_correlation_function))
    
    y_label <- 'p21-MS2 autocorrelation'
    y_axis_breaks <- c(-0.5, 0, 0.5, 1)
    y_axis_limits <- c(-0.5, 1)
    x_axis_breaks <- c(0, 5, 10, 20, 40)
    
  } else if (correlation_function_name == 'p53_p21rna_ccf') {
    
    delay_times <- seq(-(dim(sim_correlation_function)[2] - 1)/2*(p21rna_meas$time[2] - p21rna_meas$time[1]),
                       (dim(sim_correlation_function)[2] - 1)/2*(p21rna_meas$time[2] - p21rna_meas$time[1]),
                       p21rna_meas$time[2] - p21rna_meas$time[1])
    meas_correlation_function <- data.frame(t(compute_correlation_function(p53_meas, p21rna_meas)))
    colnames(meas_correlation_function) <- gsub('V', 'p53_p21rna_ccf_', colnames(meas_correlation_function))
    
    y_label <- 'p53 - p21-MS2 crosscorrelation'
    y_axis_breaks <- c(-0.2, 0, 0.2, 0.4)
    y_axis_limits <- c(-0.2, 0.4)
    x_axis_breaks <- c(-40, -20, -10, -5, 0, 5, 10, 20, 40)
  }
  
  correlation_function_plot_data <-
    data.frame(delay_times,
               c(t(meas_correlation_function)),
               c(t(sim_correlation_function[
                 which.min(rowSums(sweep(sim_correlation_function, 2, c(t(meas_correlation_function)))^2)),])))
  colnames(correlation_function_plot_data) <- 
    c('delay', 'meas_correlation_function', 'sim_correlation_function')
  
  print(which.min(rowSums(sweep(sim_correlation_function, 2, c(t(meas_correlation_function)))^2)))
  
  all_sims <- melt(t(sim_correlation_function))
  colnames(all_sims)[1:2] <- c('delay', 'cell')
  all_sims$delay <- 
    rep(delay_times, length(unique(all_sims$cell)))
  
  sim_upper_bound <- array(dim = length(unique(all_sims$delay)))
  sim_lower_bound <- array(dim = length(unique(all_sims$delay)))
  
  for (timepoint in 1:length(unique(all_sims$delay))) {
    sim_upper_bound[timepoint] <- 
      quantile(all_sims[all_sims$delay == unique(all_sims$delay)[timepoint], 'value'], 0.975)
    sim_lower_bound[timepoint] <- 
      quantile(all_sims[all_sims$delay == unique(all_sims$delay)[timepoint], 'value'], 0.025)
  }
  
  correlation_function_plot_data$sim_upper_bound <- sim_upper_bound
  correlation_function_plot_data$sim_lower_bound <- sim_lower_bound
  
  correlation_function_plot <- ggplot(correlation_function_plot_data, aes(x = delay)) +
    geom_ribbon(aes(ymin = sim_lower_bound, ymax = sim_upper_bound), fill = ribbon_color, alpha = 0.5) +
    geom_line(aes(y = meas_correlation_function, color = 'Data'), size = 1.5) +
    geom_line(aes(y = sim_correlation_function, color = model_label), size = 1.5) +
    scale_color_manual(values = c('black', line_color)) +
    xlab('Delay (h)') + ylab(y_label) +
    theme(legend.title = element_blank(), legend.justification = c(1, 1), legend.position = c(1, 1),
          #legend.box.margin = ggplot2::margin(rep(10, 4))) +
          legend.box.margin = ggplot2::margin(rep(20, 4))) +
    scale_x_continuous(trans = pseudo_log_trans(), breaks = x_axis_breaks) +
    scale_y_continuous(breaks = y_axis_breaks) +
    coord_cartesian(ylim = y_axis_limits) 
  print(correlation_function_plot)
  
  return(correlation_function_plot)
}
