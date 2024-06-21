# For running on cluster
#setwd('/michorlab/jdean/p21_central_dogma')
#args <- commandArgs(TRUE)

# For running on local computer
setwd('/Users/jamiedean/Documents/Papers/Code for Publications/p21-iffl-regulation')
args <- c('test', 0.5, 1000, NA, NA, NA, NA, 50, 5, 0.3, 0.3)

source('code/analysis_functions/preprocessing.R')
source('code/stochastic_models/p21rna_ssa.R')

################################################################################################################

model <- 'telegraph_standard'
#model <- 'telegraph_iffl'

#experiment_meas <- 'IR_2min_resolution'
#experiment_meas <- 'IR_15min_resolution'
experiment_meas <- 'IR_Nutlin'

p53 <- load_data(experiment_meas, 'traces_CFP', FALSE)

if (model == 'telegraph_standard') {
  
  # args[1] is parameter set ID number
  parameters <- unlist(data.frame(
    k_on = as.numeric(args[2]),
    k_off = as.numeric(args[3]),
    alpha_RNA = as.numeric(args[8]),
    beta_RNA = as.numeric(args[9])))
  
} else if (model == 'telegraph_iffl') {
  
  # args[1] is parameter set ID number
  parameters <- unlist(data.frame(
    k_on = as.numeric(args[2]),
    k_off = as.numeric(args[3]),
    alpha_RNA = as.numeric(args[8]),
    beta_RNA = as.numeric(args[9]),
    alpha_Repressor = as.numeric(args[10]),
    beta_Repressor = as.numeric(args[11])))
}

simulated_dataset <- run_simulation(model, p53, parameters)

write.csv(simulated_dataset, paste(paste(paste(paste(
  'outputs/simulated_datasets/', experiment_meas, sep = ''), 
  model, sep = '/'),
  args[1], sep = '_simulated_dataset_'),
  '.csv', sep =''))
