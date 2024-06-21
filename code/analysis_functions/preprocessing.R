library(reshape2)

################################################################################################################

load_data <- function(experiment, filename, normalize) {
  
  raw_data <-
    read.csv(paste(paste(paste(paste('data/microscopy/', experiment, sep = ''),
                               '/', sep = ''), filename, sep = ''), '.csv', sep = ''), header = FALSE,
             na.strings = -1)

  # p21 signal for IR + Nutlin experiment is sampled at a different frequencies so remove columns with all NaNs
  if (experiment == 'IR_Nutlin' & filename == 'traces_Texas') {
    
    raw_data <- raw_data[, seq(1, dim(raw_data)[2], 3)]
  }
  
  # Select only cells with complete trace
  complete_cases <- complete.cases(raw_data)
  raw_data <- raw_data[complete_cases,]

  if (experiment == 'IR_2min_resolution' | 
      experiment == 'Simulated_telegraph_2min' | 
      experiment == 'Simulated_transcription_telegraph_2min' | 
      experiment == 'Simulated_p53_2min') {
    
    imaging_times <- seq(2, 902, 2)/60
    
  }  else if (experiment == 'IR_15min_resolution' | 
              experiment == 'Simulated_p53_15min') {
    
    imaging_times <- seq(15, 2625, 15)/60
    annotation <-
      read.csv(paste(paste('data/microscopy/', experiment, sep = ''), '/annotation.csv', sep = ''),
               header = FALSE, na.strings = -1)
    # Select only cells with complete trace
    annotation <- annotation[complete_cases,]

    condition <- array(dim = length(annotation))
    condition[annotation == "'M20170724_4_10Gy'"] <- '10 Gy'
    condition[annotation == "'M20170724_5_5Gy'"] <- '5 Gy'
    condition[annotation == "'M20170724_6_2point5Gy'"] <- '2.5 Gy'
    
  } else if (experiment == 'IR_Nutlin') {
    
    if (filename == 'traces_Texas') {
      
      imaging_times <- seq(5, 1285, 15)/60 
      
    } else {
      
      imaging_times <- seq(5, 1285, 5)/60 
    }
    
    annotation <-
      read.csv(paste(paste('data/microscopy/', experiment, sep = ''), '/annotation.csv', sep = ''),
               header = FALSE, na.strings = -1)
    # Select only cells with complete trace
    annotation <- annotation[complete_cases,]
    
    condition <- array(dim = length(annotation))
    condition[annotation == "'M20171129_nutlin_10uM'"] <- 'IR + nutlin'
    condition[annotation == "'M20171129_nonTreated'"] <- 'IR'
    
  } else if (experiment == 'IR_p53siRNA_p21siRNA') {
    
    imaging_times <- seq(0.25, 45, 0.25)
    
    condition <-
      read.csv(paste(paste('data/microscopy/', experiment, sep = ''), '/group_number.csv',
                     sep = ''), header = FALSE, na.strings = -1)
    # Select only cells with complete trace
    condition <- condition[complete_cases,]
    
    condition[condition == 1] <- 'p53 siRNA p21 siRNA'
    condition[condition == 2] <- 'Control siRNA'
    condition[condition == 3] <- 'p21 siRNA'
    
  } 
  # For MS2 traces, perform baseline subtraction for each cell individually
  if (filename == 'traces_MS2' |
      (experiment == 'IR_Nutlin' & filename == 'traces_YFP')) {
    
    df <- data.frame(t((t(as.matrix(data.frame(t(raw_data)))) - 
                          apply(as.matrix(data.frame(t(raw_data))), 2, min))))
    
  } else {
    
    df <- data.frame(t(raw_data)) 
  }

  df$time <- imaging_times
  df_melted <- melt(df, id.vars = 'time')
  
  if (experiment == 'IR_15min_resolution' | 
      experiment == 'Simulated_damage_divisions' | 
      experiment == 'Simulated_p53_15min' |
      experiment == 'IR_Nutlin') {
    
    df_melted$condition <- rep(condition, each = dim(raw_data)[2])
    
  } else if (experiment == 'IR_p53siRNA_p21siRNA') {
    
    df_melted$condition <- rep(condition, each = dim(raw_data)[2])
  }

  print(c('Total cells:', length(unique(df_melted$variable))))

  if (filename != 'traces_area' & filename != 'traces_solidity' &
      experiment != 'Simulated_telegraph_2min' &
      experiment != 'Simulated_transcription_telegraph_2min' & 
      experiment != 'Simulated_p53_2min' & 
      experiment != 'Simulated_p53_15min') {
    
    # Remove cells that have nuclear segmentation errors (cells that do not have complete area traces)
    df_melted <- remove_cells_with_segmentation_errors(experiment, df_melted)
    
  }
  
  # Remove cells that have extremely high p53 or p21 signals
  if (filename == 'traces_CFP' | filename == 'traces_Texas') {
    
    df_melted <- remove_outlier_cells(df_melted) 
  }
  
  # Scale fluorescent signal to protein copy number (approximately) based on Burgin et al. Analyst 2014
  if ((filename == 'traces_CFP' | filename == 'traces_Texas') & experiment == 'IR_2min_resolution') {
    
    df_melted$value <- df_melted$value*30
    
  } else if ((filename == 'traces_CFP' | filename == 'traces_Texas') & experiment == 'IR_15min_resolution') {
    
    df_melted$value <- df_melted$value*5
    
  } else if ((filename == 'traces_CFP' | filename == 'traces_Texas') & experiment == 'IR_Nutlin') {
    
    df_melted$value <- df_melted$value*42
  }
  
  # Normalize signal by dividing by maximum signal in dataset (after outliers have been removed)
  if (normalize == TRUE & filename != 'traces_MS2') {
    
    df_melted$value <- df_melted$value/max(df_melted$value)
  }
  
  return(df_melted)
}

remove_outlier_cells <- function(traces) {
  # Remove cells with extremely high p53 or p21 signals
  # Outliers defined as cells with p53 signal 2 x the interquartile range (IQR) higher than the 75th quantile
  # on the log base 10 scale
  # Note that 2 x IQR is conservative (c.f. R boxplot uses 1.5 x IQR)
  
  if(any(colnames(traces) == 'condition')) {
    
    outlier_cells <- c()
    
    for (d in 1:length(unique(traces$condition))) {
      
      single_condition <- unique(traces$condition)[d]
      traces_single_condition <- traces[traces$condition == single_condition,]
      iqr <- IQR(log10(traces_single_condition$value))
      quantiles <- quantile(log10(traces_single_condition$value), probs = c(0.25, 0.75))
      print(as.character(unique(traces_single_condition[log10(traces_single_condition$value) > quantiles[2] + 2*iqr,
                                                  'variable'])))
      outlier_cells <- 
        append(outlier_cells,
               as.character(
                 unique(traces_single_condition[log10(traces_single_condition$value) > quantiles[2] + 2*iqr, 'variable'])))
    }
  } else {
    
    iqr <- IQR(log10(traces$value))
    quantiles <- quantile(log10(traces$value), probs = c(0.25, 0.75))
    outlier_cells <- unique(traces[log10(traces$value) > quantiles[2] + 2*iqr, 'variable'])
  }
  
  traces <- traces[!(traces$variable %in% outlier_cells),]
  
  print(c('Removed outlier cells:', outlier_cells))
  
  return(traces)
}

remove_cells_with_segmentation_errors <- function(experiment, traces_melted) {
  # Remove cells that do not have area or solidity data

  area <- load_data(experiment, 'traces_area', FALSE)
  solidity <- load_data(experiment, 'traces_solidity', FALSE)
  segmented_cells <- intersect(unique(area$variable), unique(solidity$variable))
  traces_melted <- traces_melted[traces_melted$variable %in% segmented_cells,]

  return(traces_melted)
}
