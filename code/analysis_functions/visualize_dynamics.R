library(circlize)
library(ComplexHeatmap)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(gridExtra)
library(reshape2)

source('code/analysis_functions/custom_plot_theme.R')

################################################################################################################

plot_single_cell_dynamics <- function(p53_traces, p21rna_traces, p21_traces, cell) {
  
  p53_data <- p53_traces[p53_traces$variable == cell,]
  
  if (is.null(p21rna_traces) == FALSE) {
    p21rna_data <- p21rna_traces[p21rna_traces$variable == cell,]
  } else {
    p21rna_data <- p53_traces[p53_traces$variable == cell,]
    p21rna_data$value <- NA
  }
  
  if (is.null(p21_traces) == FALSE) {
    p21_data <- p21_traces[p21_traces$variable == cell,]
  } else {
    p21_data <- p53_traces[p53_traces$variable == cell,]
    p21_data$value <- NA
  }
    
  p_data <- data.frame(p53_data$time, p53_data$value, p21rna_data$time, p21rna_data$value, p21_data$time,
                       p21_data$value)
  colnames(p_data) <- c('p53_time', 'p53_signal', 'p21rna_time', 'p21rna_signal', 'p21_time', 'p21_signal')
    
  p_p53 <- ggplot(p_data) + 
    geom_line(aes(x = p53_time, y = p53_signal)) +
    ylab('p53 protein (a.u.)') + xlab('Time (h)')
  p_p21rna <- ggplot(p_data) + 
    geom_line(aes(x = p21rna_time, y = p21rna_signal)) +
    ylab('p21 RNA (a.u.)') + xlab('Time (h)')
  p_p21 <- ggplot(p_data) + 
    geom_line(aes(x = p21_time, y = p21_signal)) +
    ylab('p21 protein (a.u.)') + xlab('Time (h)')
  p_combined <- ggarrange(p_p53, p_p21rna, p_p21, ncol = 1, align = 'v')
  print(p_combined)
  
}

plot_summary_dynamics <- function(traces, log_transform_expression = FALSE) {
  
  time <- unique(traces$time)
  
  mean <- array(dim = length(time))
  sd <- array(dim = length(time))
  coef_var <- array(dim = length(time))
  fano_factor <- array(dim = length(time))
  
  if (log_transform_expression == TRUE) {
    traces$value <- log(traces$value)
  }
  
  for (t in c(1:length(time))) {
    
    mean[t] <- mean((traces$value[traces$time == time[t]]), na.rm = TRUE)
    sd[t] <- sd((traces$value[traces$time == time[t]]), na.rm = TRUE)
    coef_var[t] <- sd[t]/mean[t]
    fano_factor[t] <- sd[t]^2/mean[t]
  }
  
  summary_dynamics <- data.frame(time, mean, sd, coef_var, fano_factor)
  
  p_mean <- ggplot(summary_dynamics, aes(time)) + 
    geom_line(aes(y = mean), size = 2) + 
    ylab('Mean') + xlab('time (h)') +
    custom_plot_theme

  p_sd <- ggplot(summary_dynamics, aes(time)) + 
    geom_line(aes(y = sd), size = 2) + 
    ylab('Standard deviation') + xlab('time (h)') +
    custom_plot_theme
  
  p_coef_var <- ggplot(summary_dynamics, aes(time)) + 
    geom_line(aes(y = coef_var), size = 2) + 
    ylab('Coefficient of variation') + xlab('time (h)') +
    custom_plot_theme
  
  p_fano_factor <- ggplot(summary_dynamics, aes(time)) + 
    geom_line(aes(y = fano_factor), size = 2) + 
    ylab('Fano factor') + xlab('time (h)') +
    custom_plot_theme
  
  p_summary_dynamics <- grid.arrange(p_mean, p_sd, p_coef_var, p_fano_factor)
  print(p_summary_dynamics)
  
  return(summary_dynamics)
}

plot_summary_dynamics_by_condition <- function(traces, experiment, log_transform_expression = FALSE, return_plot = FALSE) {
  
  if (experiment == 'IR_15min_resolution') {
    
    summary_dynamics_10Gy <- plot_summary_dynamics(traces[traces$condition == '10 Gy',], log_transform_expression)
    summary_dynamics_5Gy <- plot_summary_dynamics(traces[traces$condition == '5 Gy',], log_transform_expression)
    summary_dynamics_2.5Gy <- plot_summary_dynamics(traces[traces$condition == '2.5 Gy',], log_transform_expression)
    
    summary_dynamics_10Gy$condition <- '10 Gy'
    summary_dynamics_5Gy$condition <- '5 Gy'
    summary_dynamics_2.5Gy$condition <- '2.5 Gy'
    
    summary_dynamics <- rbind(summary_dynamics_10Gy, summary_dynamics_5Gy, summary_dynamics_2.5Gy)
    summary_dynamics$condition <- factor(summary_dynamics$condition, levels = c('2.5 Gy', '5 Gy', '10 Gy'))
    
  } else if (experiment == 'IR_Nutlin') {
    
    summary_dynamics_0uM <- plot_summary_dynamics(traces[traces$condition == 'IR',], log_transform_expression)
    summary_dynamics_10uM <- plot_summary_dynamics(traces[traces$condition == 'IR + nutlin',], log_transform_expression)
    
    summary_dynamics_0uM$condition <- 'IR'
    summary_dynamics_10uM$condition <- 'IR + nutlin'
    
    summary_dynamics <- rbind(summary_dynamics_0uM, summary_dynamics_10uM)
    summary_dynamics$condition <- factor(summary_dynamics$condition, levels = c('IR', 'IR + nutlin'))
  }
  
  p_mean <- ggplot(summary_dynamics) +
    geom_line(aes(x = time, y = mean, colour = condition), size = 2) +
    geom_ribbon(aes(x = time, ymin = mean - sd, ymax = mean + sd, fill = condition), alpha = 0.5,
                show.legend = FALSE) +
    xlab('Time (h)') +
    ylab('Mean (a.u.)') +
    #ylab('Nuclear area (a.u.)') +
    custom_plot_theme + 
    theme(legend.position = c(0.25, 0.87))

  p_sd <- ggplot(summary_dynamics) + 
    geom_line(aes(x = time, y = sd, colour = condition), size = 2) + 
    ylab('Standard deviation (a.u.)') + xlab('Time (h)') +
    custom_plot_theme
  
  p_coef_var <- ggplot(summary_dynamics) + 
    geom_line(aes(x = time, y = coef_var, colour = condition), size = 2) + 
    ylab('Coefficient of variation') + xlab('Time (h)') +
    custom_plot_theme
  
  p_fano_factor <- ggplot(summary_dynamics) + 
    geom_line(aes(x = time, y = fano_factor, colour = condition), size = 2) + 
    ylab('Fano factor') + xlab('Time (h)') +
    custom_plot_theme
  
  print(p_mean)
  print(p_sd)
  print(p_coef_var)
  print(p_fano_factor)

  if (return_plot == TRUE) {
    
    return(p_mean)
    
  } else {
   
    return(summary_dynamics)
  }
}

plot_single_cell_dynamics_heatmaps <- function(traces, signal_name, column_split = FALSE) {
  
  heatmaps <- list()
  
  if (signal_name == 'p53') {
    plot_color <- plot_colors[1]
  } else if (signal_name == 'p21-MS2') {
    plot_color <- plot_colors[2]
  } else if (signal_name == 'p21') {
    plot_color <- plot_colors[3]
  }
  
  traces_wide <- acast(traces, variable ~ time, value.var = 'value')
  
  if (min(traces_wide) <= 0) {
    traces_wide <- traces_wide + min(traces_wide) + 1
  }
  
  traces_wide <- log2(traces_wide/traces_wide[, 1])
  quantile_5 <- quantile(traces_wide, 0.05, na.rm = TRUE)
  quantile_95 <- quantile(traces_wide, 0.95, na.rm = TRUE)
  
  for (single_condition in 1:length(unique(traces$condition))) {
    
    traces_wide <- 
      acast(traces[traces$condition == unique(traces$condition)[single_condition],], variable ~ time, value.var = 'value')

    if (min(traces_wide) <= 0) {
      traces_wide <- traces_wide + min(traces_wide) + 1
    }
    
    traces_wide <- log2(traces_wide/traces_wide[, 1])
    
    x_axis_labels <- rep('', dim(traces_wide)[2])
    x_axis_labels_display <- 
      which(sapply((unique(traces$time) - traces$time[1])/3, function(i) i == as.integer(i)))
    x_axis_labels[x_axis_labels_display] <- unique(traces$time)[x_axis_labels_display] - traces$time[1]
    
    
    if (column_split == TRUE) {
      
      sirna_administration_frame <- 80
      column_split_custom <- 
        data.frame(factor(rep(c('Pre-siRNA', 'Post-siRNA'),
                              c(sirna_administration_frame, ncol(traces_wide) - sirna_administration_frame)),
                          levels = c('Pre-siRNA', 'Post-siRNA')))  
      
    } else if (column_split == FALSE) {
      
      column_split_custom <- NULL
    }
    
    heatmaps[single_condition] <-
      Heatmap(
        traces_wide,
        col = colorRamp2(c(quantile_5, quantile_95), c('white', plot_color)),
        column_order = seq(1, dim(traces_wide)[2]),
        row_order = order(rowSums(traces_wide)),
        show_row_names = FALSE,
        column_labels = x_axis_labels,
        show_column_names = TRUE,
        column_names_rot = 0,
        column_names_side = 'top',
        border = TRUE,
        name = signal_name,
        column_title = 'Time (h)',
        column_split = column_split_custom,
        row_title = paste(unique(traces$condition)[single_condition], '\n \n Cell'),
        column_title_gp = grid::gpar(fontsize = 24),
        row_title_gp = grid::gpar(fontsize = 24),
        column_names_gp = grid::gpar(fontsize = 24),
        row_names_gp = grid::gpar(fontsize = 24),
        heatmap_legend_param = list(
          title_gp = gpar(fontsize = 24),
          labels_gp = gpar(fontsize = 24))
        )
  }
  
  if (length(unique(traces$condition)) == 2) {
    
    heatmaps <- draw(heatmaps[[1]] %v% heatmaps[[2]], ht_gap = unit(1, 'cm'))
    
  } else if (length(unique(traces$condition)) == 3) {
    
    heatmaps <- draw(heatmaps[[1]] %v% heatmaps[[2]] %v% heatmaps[[3]], ht_gap = unit(1, 'cm'))
  }
  
  return(heatmaps)
}

plot_change_in_traces <- function(traces) {
  
  traces_wide_0uM <- acast(traces[traces$condition == 'IR',], variable ~ time, value.var = 'value')
  traces_wide_10uM <- acast(traces[traces$condition == 'IR + nutlin',], variable ~ time, value.var = 'value')
  
  plot_data <- 
    data.frame(
      traces = c(traces_wide_0uM[, 1], traces_wide_0uM[, dim(traces_wide_0uM)[2] - 4], 
                 traces_wide_10uM[, 1], traces_wide_10uM[, dim(traces_wide_10uM)[2] - 4]),
      time = c(rep('0 h', dim(traces_wide_0uM)[1]), rep('21 h', dim(traces_wide_0uM)[1]),
               rep('0 h', dim(traces_wide_10uM)[1]), rep('21 h', dim(traces_wide_10uM)[1])),
      condition = c(rep('IR', 2*dim(traces_wide_0uM)[1]), rep('IR + nutlin', 2*dim(traces_wide_10uM)[1])),
      cell = c(rep(rownames(traces_wide_0uM), 2), rep(rownames(traces_wide_10uM), 2)))
  
  p_change <- ggpaired(plot_data, x = 'time', y = 'traces', id = 'cell',
                       fill = 'condition',
                       facet.by = 'condition',
                       line.color = 'gray', line.size = 0.05,
                       legend = 'none',
                       xlab = 'Time', ylab = 'Nuclear area (a.u.)',
                       point.size = 0.5,
                       ggtheme = custom_plot_theme) + 
    stat_compare_means(paired = TRUE, method = 't.test', label = 'p.format', hjust = -0.15) +
    scale_color_manual(values = c(plot_colors[1], plot_colors[2]))
  
  return(p_change)  
}
