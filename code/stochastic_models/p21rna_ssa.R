# Note that ssar library is the version at https://github.com/jamiedean/ssar
#devtools::install_github('jamiedean/ssar')

library(ggplot2)
library(ssar)

# Continuous time discrete state Markov model with time-dependent propensities

################################################################################################################

run_simulation_p21rna_telegraph_standard_single_cell <- function(p53_data, params, cell_number) {
  # Random telegraph model
  
  p53_data_shifted <- p53_data[p53_data$time > 50/60,]
  
  p53_cell <- p53_data_shifted[p53_data_shifted$variable == unique(p53_data_shifted$variable)[cell_number],]
  
  # Initial conditions
  X <- matrix(c(Gene_off = 1, Gene_on = 0, RNA = 0), ncol = 3)
  
  # Reaction matrix (rows - states; columns - reactions)
  v <- matrix(c(+1, -1, 0,  0,#   Promoter off state
                -1, +1, 0,  0,#   Promoter on state
                0,  0,  +1, -1),# Nascent RNA
              nrow = 3, byrow = TRUE)
  
  # Propensity vector
  pfun <- function(t, X, params) {
    cbind(
      params['k_off']*X[, 2],
      params['k_on']*p53_cell[which.min(abs(p53_cell$time - t)), 'value']*X[, 1],
      params['alpha_RNA']*X[, 2],
      params['beta_RNA']*X[, 3]
    )
  }
  simulation <- 
    ssa(X, pfun, v, params, tmin = min(p53_cell$time), tmax = max(p53_cell$time), nsim = 1, plot.sim = FALSE)
  
  observed_sim_data <- sapply(1:length(p53_cell$time),
                              function(x) simulation$Var3[which.min(abs(simulation$Time - p53_cell$time[x]))])
  
  sim_output <- data.frame(p53_cell$time, p53_cell$variable, observed_sim_data)
  colnames(sim_output) <- c('time', 'variable', 'value')
  
  return(sim_output)
}

run_simulation_p21rna_telegraph_iffl_single_cell <- function(p53_data, params, cell_number) {
  # Random telegraph model with incoherent feedfoward loop
  
  p53_data_shifted <- p53_data[p53_data$time > 50/60,]
  
  p53_cell <- p53_data_shifted[p53_data_shifted$variable == unique(p53_data_shifted$variable)[cell_number],]
  
  repressor_initial <- p53_data[p53_data$variable == unique(p53_data$variable)[cell_number], 'value'][1]
  
  # Initial conditions
  X <- matrix(c(Gene_off = 1, Gene_on = 0, RNA = 0, Repressor = repressor_initial), ncol = 4)
  
  # Reaction matrix (rows - states; columns - reactions)
  v <- matrix(c(+1, -1, 0,  0,  0,  0,#   Promoter off state
                -1, +1, 0,  0,  0,  0,#   Promoter on state
                0,  0,  +1, -1, 0,  0,#   Nascent RNA
                0,  0,  0,  0,  +1, -1),# Repressor
              nrow = 4, byrow = TRUE)
  
  # Propensity vector
  pfun <- function(t, X, params) {
    cbind(
      params['k_off']*X[, 2],
      params['k_on']*(
        ifelse((p53_cell[which.min(abs(p53_cell$time - t)), 'value'] - X[, 4]) > 0,
               p53_cell[which.min(abs(p53_cell$time - t)), 'value'] - X[, 4], 0))*X[, 1],
      params['alpha_RNA']*X[, 2],
      params['beta_RNA']*X[, 3],
      params['alpha_Repressor']*p53_cell[which.min(abs(p53_cell$time - t)), 'value'],
      params['beta_Repressor']*X[, 4]
    )
  }
  simulation <- 
    ssa(X, pfun, v, params, tmin = min(p53_cell$time), tmax = max(p53_cell$time), nsim = 1, plot.sim = FALSE)
  
  observed_sim_data <- sapply(1:length(p53_cell$time),
                              function(x) simulation$Var3[which.min(abs(simulation$Time - p53_cell$time[x]))])
  observed_sim_repressor <- 
    sapply(1:length(p53_cell$time), 
           function(x) simulation$Var4[which.min(abs(simulation$Time - p53_cell$time[x]))])
  
  sim_output <- data.frame(p53_cell$time, p53_cell$variable, observed_sim_data, observed_sim_repressor)
  colnames(sim_output) <- c('time', 'variable', 'value', 'repressor')
  
  return(sim_output)
}

run_simulation <- function(model, p53_data, parameters) {

  if (model == 'telegraph_standard') {
    
    simulated_data <- 
      lapply(1:length(unique(p53_data$variable)),
             function(cell_number) run_simulation_p21rna_telegraph_standard_single_cell(
               p53_data, parameters, cell_number))
    
  } else if (model == 'telegraph_iffl') {
    
    simulated_data <- 
      lapply(1:length(unique(p53_data$variable)),
             function(cell_number) run_simulation_p21rna_telegraph_iffl_single_cell(
               p53_data, parameters, cell_number))
  }
  
  simulated_data <- do.call(rbind, simulated_data)

  return(simulated_data)
}

plot_single_cell_simulation <- function(p53_data, p21rna_data, model, params, cell_number) {
  
  if (model == 'telegraph_standard') {
    
    observed_sim_data <- run_simulation_p21rna_telegraph_standard_single_cell(p53_data, params, cell_number)
    model_color <- plot_colors[1]
    
  } else if (model == 'telegraph_iffl') {
    
    observed_sim_data <- run_simulation_p21rna_telegraph_iffl_single_cell(p53_data, params, cell_number)
    model_color <- plot_colors[2]
  }
  
  p53_cell <- p53_data[p53_data$variable == unique(p53_data$variable)[cell_number],]
  p21rna_cell <- p21rna_data[p21rna_data$variable == as.character(unique(p53_data$variable)[cell_number]),]
  
  # Fluroescence maturation correction
  # CFP has ~ 50 min maturation time
  p53_cell <- p53_cell[p53_cell$time > 50/60,]
  p21rna_cell <- p21rna_cell[1:length(p53_cell$time),]
  
  data <- p53_cell
  data$p21rna_meas <- p21rna_cell$value
  data$p21rna_sim <- observed_sim_data$value
  
  meas_plot <- ggplot(data = data, aes(x = time)) +
    geom_line(aes(y = value/750), color = 'gray', size = 1) +
    geom_line(aes(y = p21rna_meas), color = 'black', size = 1) +
    xlab('Time (h)') + ylab('Fluorescence (a.u.)') +
    theme(legend.title = element_blank())
  print(meas_plot)
  
  sim_plot <- ggplot(data = data, aes(x = time)) +
    geom_line(aes(y = value/750), color = 'gray', size = 1) +
    geom_line(aes(y = p21rna_meas), color = 'black', size = 1) +
    geom_line(aes(y = p21rna_sim), color = model_color, size = 1) +
    xlab('Time (h)') + ylab('Fluorescence (a.u.)') +
    theme(legend.position = 'none')
  
  if (model == 'telegraph_iffl') {
    
    data$repressor_sim <- observed_sim_data$repressor
    sim_plot <- sim_plot + 
      #geom_line(aes(y = data$repressor_sim/750, color = 'Repressor'), size = 1) +
      scale_color_manual(labels = c('p53', 'p21 RNA simulated', 'p21 RNA measured', 'Repressor'),
                         values = c('gray', plot_colors[2], 'black', plot_colors[3]))
    
  } else {
    
    sim_plot <- sim_plot + 
      scale_color_manual(labels = c('p53', 'p21 RNA simulated', 'p21 RNA measured'),
                         values = c('gray', plot_colors[2], 'black'))
  }
  
  print(sim_plot)
  
  acf_data <- c()
  p21rna_meas_acf <- acf(data$p21rna_meas, lag.max = length(data$time), plot = FALSE)
  acf_data$delay <- c(p21rna_meas_acf$lag*(data$time[2] - data$time[1]))
  acf_data$p21rna_meas_acf <- c(p21rna_meas_acf$acf)
  acf_data$p21rna_sim_acf <- c(acf(data$p21rna_sim, lag.max = length(data$time), plot = FALSE)$acf)
  acf_data <- data.frame(acf_data)
  
  acf_plot <- ggplot(data = acf_data, aes(x = delay)) +
    geom_line(aes(y = p21rna_meas_acf), color = 'black', size = 1.5) +
    geom_line(aes(y = p21rna_sim_acf), color = model_color, size = 1.5) +
    xlab('Delay (h)') + ylab('p21-MS2 autocorrelation') +
    scale_x_continuous(trans = pseudo_log_trans(), breaks = c(0, 5, 10, 20, 40)) +
    scale_color_manual(labels = c('Experiment', 'Simulation'),
                       values = c('black', model_color)) +
    theme(legend.title = element_blank(), legend.justification = c(1, 1), legend.position = c(1, 1),
          legend.box.margin = ggplot2::margin(rep(10, 4)))
  print(acf_plot)
  
  ccf_data <- c()
  p53_p21rna_meas_ccf <- ccf(data$value, data$p21rna_meas, lag.max = length(data$time), plot = FALSE)
  ccf_data$delay <- c(p53_p21rna_meas_ccf$lag*(data$time[2] - data$time[1]))
  ccf_data$p53_p21rna_meas_ccf <- c(p53_p21rna_meas_ccf$acf)
  ccf_data$p53_p21rna_sim_ccf <- c(ccf(data$value, data$p21rna_sim, lag.max = length(data$time), plot = FALSE)$acf)
  ccf_data <- data.frame(ccf_data)
  
  ccf_plot <- ggplot(data = ccf_data, aes(x = delay)) +
    geom_line(aes(y = p53_p21rna_meas_ccf), color = 'black', size = 1.5) +
    geom_line(aes(y = p53_p21rna_sim_ccf), color = model_color, size = 1.5) +
    scale_color_manual(labels = c('Experiment', 'Simulation'),
                       values = c('black', model_color)) +
    xlab('Delay (h)') + ylab('p53 - p21-MS2 crosscorrelation') +
    scale_x_continuous(trans = pseudo_log_trans(), breaks = c(-40, -20, -10, -5, 0, 5, 10, 20, 40)) +
    theme(legend.title = element_blank(), legend.justification = c(1, 1), legend.position = c(1, 1),
          legend.box.margin = ggplot2::margin(rep(10, 4)))
  print(ccf_plot)
  
  return(list(sim_plot, acf_plot, ccf_plot))
}
