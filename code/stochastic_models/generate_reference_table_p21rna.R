# For running on cluster
#args <- commandArgs(TRUE)

# For running on local computer
args <- c('test')

################################################################################################################

setwd('/Users/jamiedean/Documents/Papers/Code for Publications/p21-iffl-regulation')
#setwd('/michorlab/jdean/p21_central_dogma')

source('code/analysis_functions/approximate_bayesian_computation.R')

################################################################################################################

model <- 'telegraph_standard'
#model <- 'telegraph_iffl'

experiment_meas <- 'IR_2min_resolution'
#experiment_meas <- 'IR_15min_resolution'
#experiment_meas <- 'IR_Nutlin'

p53 <- load_data(experiment_meas, 'traces_CFP', FALSE)

# Fluroescence maturation correction
# CFP has ~ 50 min maturation time
p53 <- p53[p53$time > 50/60,]

simulated_dataset <- 
  read.csv(
    paste(paste(paste(paste(
          'outputs/simulated_datasets/', experiment_meas, sep = ''), 
          model, sep = '/'),
          args, sep = '_simulated_dataset_'),
          '.csv', sep =''))

for (condition in unique(p53$condition)) {
  
  correlation_functions <- 
    generate_correlation_functions(
      p53[p53$condition == condition,], simulated_dataset[which(p53$condition == condition),], NULL)

  p21rna_autocorrelation <- correlation_functions[[2]]
  p53_p21rna_crosscorrelation <- correlation_functions[[4]]

  write.csv(p21rna_autocorrelation, paste(paste(paste(paste(paste(paste(
    'outputs/correlation_functions/', experiment_meas, sep = ''),
    condition, sep = '/'),
    model, sep = '/'),
    'p21rna_autocorrelation', sep = '/'),
    args, sep = '_'),
    '.csv', sep = ''), row.names = FALSE)
  write.csv(p53_p21rna_crosscorrelation, paste(paste(paste(paste(paste(paste(
    'outputs/correlation_functions/', experiment_meas, sep = ''),
    condition, sep = '/'),
    model, sep = '/'),
    'p53_p21rna_crosscorrelation', sep = '/'),
    args, sep = '_'),
    '.csv', sep = ''), row.names = FALSE)
}
