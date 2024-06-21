library(AnnotationDbi)
library(ggsurvfit)
library(org.Hs.eg.db)
library(pheatmap)
library(RadioGx)
library(readr)
library(survival)
library(survminer)
library(tidyr)

source('code/analysis_functions/custom_plot_theme.R')
source('code/analysis_functions/functional_analysis.R')

###############################################################################################################

signature_outcome_association_cell_lines_ccle <- function() {
  # Tissue-wise correlation between IR + Nutlin vs IR signature and radiosensitivity
  
  Cleveland <- downloadRSet('Cleveland', saveDir = 'data')
  
  cell_line_metadata <- cellInfo(Cleveland)
  
  auc <- sensitivityProfiles(Cleveland)['AUC_recomputed']
  auc <- auc[-which(duplicated(sensitivityInfo(Cleveland)[, 'sampleid'])),]
  names(auc) <- sensitivityInfo(Cleveland)[-which(duplicated(sensitivityInfo(Cleveland)[, 'sampleid'])), 'sampleid']
  
  mprofSummary <- 
    summarizeMolecularProfiles(Cleveland, mDataType = 'rnaseq', summary.stat = 'median', fill.missing = FALSE)
  
  gsva_results <- 
    compute_geneset_variation_analysis_scores(mprofSummary, 'custom', 'ir_nutlin_vs_ir_signature', FALSE, TRUE)
  
  auc_matched_rna <- auc[names(auc) %in% colnames(gsva_results)]
  
  gsva_results_transposed <- t(gsva_results)
  colnames(gsva_results_transposed) <- c('ir_nutlin_vs_ir_up', 'ir_nutlin_vs_ir_down')
  cell_line_metadata_gsva <- merge(cell_line_metadata, gsva_results_transposed, by = 'row.names')
  row.names(cell_line_metadata_gsva) <- cell_line_metadata_gsva$Row.names
  cell_line_metadata_gsva_auc <- merge(cell_line_metadata_gsva, data.frame(auc_matched_rna), by = 'row.names')
  
  signature_auc_lm_coefficient <- matrix(nrow = length(unique(cell_line_metadata_gsva_auc$Primarysite)), ncol = 5)
  rownames(signature_auc_lm_coefficient) <- unique(cell_line_metadata_gsva_auc$Primarysite)
  
  for (i in 1:length(unique(cell_line_metadata_gsva_auc$Primarysite))) {
    
    primary_site <- unique(cell_line_metadata_gsva_auc$Primarysite)[i]
    cell_line_metadata_gsva_auc_primary_site <-
      cell_line_metadata_gsva_auc[cell_line_metadata_gsva_auc$Primarysite == primary_site,]
    
    if (nrow(cell_line_metadata_gsva_auc_primary_site) >= 10) {
      
      lm_ir_nutlin_vs_ir_up <- 
        summary(lm(auc_matched_rna ~ ir_nutlin_vs_ir_up,
                   data = cell_line_metadata_gsva_auc_primary_site))['coefficients']$coefficients
      lm_ir_nutlin_vs_ir_down <- 
        summary(lm(auc_matched_rna ~ ir_nutlin_vs_ir_down, 
                   data = cell_line_metadata_gsva_auc_primary_site))['coefficients']$coefficients
      
      signature_auc_lm_coefficient[i,] <- 
        c(lm_ir_nutlin_vs_ir_up['ir_nutlin_vs_ir_up', 'Estimate'],
          lm_ir_nutlin_vs_ir_up['ir_nutlin_vs_ir_up', 'Std. Error'],
          lm_ir_nutlin_vs_ir_down['ir_nutlin_vs_ir_down', 'Estimate'],
          lm_ir_nutlin_vs_ir_down['ir_nutlin_vs_ir_down', 'Std. Error'],
          nrow(cell_line_metadata_gsva_auc_primary_site))
      
    } else {
      
      signature_auc_lm_coefficient[i,] <- rep(NA, 5)
    }
  }
  
  signature_auc_lm_coefficient <- data.frame(signature_auc_lm_coefficient)
  signature_auc_lm_coefficient <- signature_auc_lm_coefficient[complete.cases(signature_auc_lm_coefficient),]
  colnames(signature_auc_lm_coefficient) <- 
    c('up_estimate', 'up_std_error', 'down_estimate', 'down_std_error', 'n_cell_lines')
  
  # Inverse-variance weighted mean
  signature_auc_lm_coefficient$up_inverse_variance_weighting_numerator <- 
    signature_auc_lm_coefficient$up_estimate*(1/signature_auc_lm_coefficient$up_std_error^2)
  signature_auc_lm_coefficient$down_inverse_variance_weighting_numerator <- 
    signature_auc_lm_coefficient$down_estimate*(1/signature_auc_lm_coefficient$down_std_error^2)
  signature_auc_lm_coefficient$up_inverse_variance_weighting_denominator <- 
    1/signature_auc_lm_coefficient$up_std_error^2
  signature_auc_lm_coefficient$down_inverse_variance_weighting_denominator <- 
    1/signature_auc_lm_coefficient$down_std_error^2 
  
  inverse_variance_weighted_mean_up <- 
    sum(signature_auc_lm_coefficient$up_inverse_variance_weighting_numerator)/
    sum(signature_auc_lm_coefficient$up_inverse_variance_weighting_denominator)
  inverse_variance_weighted_mean_down <-
    sum(signature_auc_lm_coefficient$down_inverse_variance_weighting_numerator)/
    sum(signature_auc_lm_coefficient$down_inverse_variance_weighting_denominator)
  # Standard error of inverse-variance weighted mean
  inverse_variance_weighted_std_error_up <- 
    sqrt(1/sum(signature_auc_lm_coefficient$up_inverse_variance_weighting_denominator))
  inverse_variance_weighted_std_error_down <- 
    sqrt(1/sum(signature_auc_lm_coefficient$down_inverse_variance_weighting_denominator))
  
  signature_auc_lm_coefficient['Weighted Mean',] <- 
    c(inverse_variance_weighted_mean_up, inverse_variance_weighted_std_error_up,
      inverse_variance_weighted_mean_down, inverse_variance_weighted_std_error_down,
      sum(signature_auc_lm_coefficient$n_cell_lines), rep(NA, 4))
  
  plot_data <- 
    data.frame(rep(rownames(signature_auc_lm_coefficient), 2),
               rep(signature_auc_lm_coefficient$n_cell_lines, 2),
               c(rep('Upregulated Signature', nrow(signature_auc_lm_coefficient)),
                 rep('Downregulated Signature', nrow(signature_auc_lm_coefficient))),
               c(signature_auc_lm_coefficient$up_estimate, signature_auc_lm_coefficient$down_estimate),
               c(signature_auc_lm_coefficient$up_std_error, signature_auc_lm_coefficient$down_std_error))
  plot_data$shape <- 16
  colnames(plot_data) <- c('primary_site', 'n_cell_lines', 'signature', 'coefficient', 'std_error', 'shape')
  plot_data[plot_data$primary_site == 'Weighted Mean', 'shape'] <- 18
  
  p <- ggplot(plot_data, aes(x = coefficient, y = factor(primary_site, ordered = TRUE), color = signature)) +
    geom_vline(xintercept = 0, color = 'gray') +
    geom_point(aes(size = n_cell_lines, shape = shape)) + scale_shape_identity() +
    geom_errorbar(aes(xmin = coefficient - std_error, xmax = coefficient + std_error), linewidth = 1) +
    facet_wrap(~ signature) +
    xlab('Coefficient') + ylab('Primary Site') +
    custom_plot_theme +
    guides(size = guide_legend('Number of cell lines'), color = 'none')
  
  z_statistic_upregulated_signature <- 
    plot_data[plot_data$primary_site == 'Weighted Mean' & plot_data$signature == 'Upregulated Signature',
              'coefficient']/
    plot_data[plot_data$primary_site == 'Weighted Mean' & plot_data$signature == 'Upregulated Signature',
              'std_error']
  p_value_upregulated_signature <- 2*pnorm(q = abs(z_statistic_upregulated_signature), lower.tail = FALSE)
  print(p_value_upregulated_signature)
    
  z_statistic_downregulated_signature <- 
    plot_data[plot_data$primary_site == 'Weighted Mean' & plot_data$signature == 'Downregulated Signature',
              'coefficient']/
    plot_data[plot_data$primary_site == 'Weighted Mean' & plot_data$signature == 'Downregulated Signature',
              'std_error']
  p_value_downregulated_signature <- 2*pnorm(q = abs(z_statistic_downregulated_signature), lower.tail = FALSE)
  print(p_value_downregulated_signature)
  
  return(p)
}

signature_outcome_association_clinical_metabric <- function(gsva, treatment = 'chemo_radiotherapy') {
  # Note: have to use do.call(survfit, list(...)) to run survfit from within a function
  
  # Clinical data
  metabric_clinical <- read_tsv('data/radiosensitivity/brca_metabric/data_clinical_patient.txt', comment = '#')
  metabric_sample <- read_tsv('data/radiosensitivity/brca_metabric/data_clinical_sample.txt', comment = '#')
  metabric_clinical <- inner_join(metabric_clinical, metabric_sample, by = ('PATIENT_ID'))
  
  # Combine data
  metabric_clinical_gsva <- 
    inner_join(metabric_clinical, gsva,  by = c('PATIENT_ID' = 'Sample_ID'))
  
  if (treatment == 'no_adjuvant_therapy') {
    
    radiotherapy <- 'NO'
    chemotherapy <- 'NO'
    hormone_therapy <- 'NO'
    
  } else if (treatment == 'chemo_radiotherapy') {
    
    radiotherapy <- 'YES'
    chemotherapy <- 'YES'
    hormone_therapy <- 'NO'
  }
  
  # Upregulated Signature
  # Cox model
  cox_model_rfs_up <- 
    coxph(Surv(RFS_MONTHS, RFS_STATUS == '1:Recurred') ~ ir_nutlin_vs_ir_up, 
          data = metabric_clinical_gsva %>%
            filter(RADIO_THERAPY == radiotherapy, CHEMOTHERAPY == chemotherapy, HORMONE_THERAPY == hormone_therapy))
  cox_model_pvalue_rfs_up <- summary(cox_model_rfs_up)$coefficients[5]
  
  cox_model_os_up <- 
    coxph(Surv(OS_MONTHS, OS_STATUS == '1:DECEASED') ~ ir_nutlin_vs_ir_up, 
          data = metabric_clinical_gsva %>%
            filter(RADIO_THERAPY == radiotherapy, CHEMOTHERAPY == chemotherapy, HORMONE_THERAPY == hormone_therapy))
  cox_model_pvalue_os_up <- summary(cox_model_os_up)$coefficients[5]
  
  # Recurrence-free survival
  p_rfs_up <- ggsurvplot(
    do.call(survfit, list(Surv(RFS_MONTHS, RFS_STATUS == '1:Recurred') ~ 
              ir_nutlin_vs_ir_up > median(ir_nutlin_vs_ir_up),
            data = metabric_clinical_gsva %>% 
              filter(RADIO_THERAPY == radiotherapy, CHEMOTHERAPY == chemotherapy, HORMONE_THERAPY == hormone_therapy))),
    pval = signif(cox_model_pvalue_rfs_up, 2),
    size = 2) +
    xlab('Time (months)') + ylab('Recurrence-Free Survival') +
    custom_plot_theme
  
  # Overall survival
  p_os_up <- ggsurvplot(
    do.call(survfit, list(Surv(OS_MONTHS, OS_STATUS == '1:DECEASED') ~ 
              ir_nutlin_vs_ir_up > median(ir_nutlin_vs_ir_up),
            data = metabric_clinical_gsva %>% 
              filter(RADIO_THERAPY == radiotherapy, CHEMOTHERAPY == chemotherapy, HORMONE_THERAPY == hormone_therapy))),
    pval = signif(cox_model_pvalue_os_up, 2),
    size = 2) +
    xlab('Time (months)') + ylab('Overall Survival') +
    custom_plot_theme
  
  # Upregulated Signature
  # Cox model
  cox_model_rfs_down <- 
    coxph(Surv(RFS_MONTHS, RFS_STATUS == '1:Recurred') ~ ir_nutlin_vs_ir_down, 
          data = metabric_clinical_gsva %>%
            filter(RADIO_THERAPY == radiotherapy, CHEMOTHERAPY == chemotherapy, HORMONE_THERAPY == hormone_therapy))
  cox_model_pvalue_rfs_down <- summary(cox_model_rfs_down)$coefficients[5]
  
  cox_model_os_down <- 
    coxph(Surv(OS_MONTHS, OS_STATUS == '1:DECEASED') ~ ir_nutlin_vs_ir_down, 
          data = metabric_clinical_gsva %>%
            filter(RADIO_THERAPY == radiotherapy, CHEMOTHERAPY == chemotherapy, HORMONE_THERAPY == hormone_therapy))
  cox_model_pvalue_os_down <- summary(cox_model_os_down)$coefficients[5]
  
  # Recurrence-free survival
  p_rfs_down <- ggsurvplot(
    do.call(survfit, list(Surv(RFS_MONTHS, RFS_STATUS == '1:Recurred') ~ 
              ir_nutlin_vs_ir_down > median(ir_nutlin_vs_ir_down),
            data = metabric_clinical_gsva %>% 
              filter(RADIO_THERAPY == radiotherapy, CHEMOTHERAPY == chemotherapy, HORMONE_THERAPY == hormone_therapy))),
    pval = signif(cox_model_pvalue_rfs_down, 2),
    size = 2) +
    xlab('Time (months)') + ylab('Recurrence-Free Survival') +
    custom_plot_theme
  
  # Overall survival
  p_os_down <- ggsurvplot(
    do.call(survfit, list(Surv(OS_MONTHS, OS_STATUS == '1:DECEASED') ~ 
              ir_nutlin_vs_ir_down > median(ir_nutlin_vs_ir_down),
            data = metabric_clinical_gsva %>% 
              filter(RADIO_THERAPY == radiotherapy, CHEMOTHERAPY == chemotherapy, HORMONE_THERAPY == hormone_therapy))),
    pval = signif(cox_model_pvalue_os_down, 2),
    size = 2) +
    xlab('Time (months)') + ylab('Overall Survival') +
    custom_plot_theme
  
  return(list(p_rfs_up, p_os_up, p_rfs_down, p_os_down))
}
