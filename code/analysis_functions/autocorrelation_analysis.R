library(gridExtra)

source('code/analysis_functions/custom_plot_theme.R')

################################################################################################################

compute_correlation_function <- function(trace1, trace2) {
  
  if (identical(trace1, trace2) == TRUE) {
    
    correlation_function <- 
      lapply(1:length(unique(trace1$variable)), function(cell_number) {
        acf(trace1[trace1$variable == unique(trace1$variable)[cell_number], 'value'],
            lag.max = length(trace1[trace1$variable == unique(trace1$variable)[cell_number], 'value']),
            plot = FALSE)$acf})
    
  } else if (identical(trace1, trace2) == FALSE) {
    
    correlation_function <- 
      lapply(1:length(unique(trace1$variable)), function(cell_number) {
        ccf(trace1[trace1$variable == unique(trace1$variable)[cell_number], 'value'],
            trace2[trace2$variable == unique(trace2$variable)[cell_number], 'value'],
            lag.max = length(trace1[trace1$variable == unique(trace1$variable)[cell_number], 'value']),
            plot = FALSE)$acf})
  }
  
  # Note that cells with constant value signal (e.g. 0) have correlation functions of NA
  # These are removed when performing the column sum but counted when dividing by the number of cells
  correlation_function <- as.data.frame(do.call(rbind, correlation_function))
  correlation_function_mean <- colSums(correlation_function, na.rm = TRUE)/nrow(correlation_function)
  
  return(correlation_function_mean)
}
