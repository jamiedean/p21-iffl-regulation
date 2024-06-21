library(AnnotationDbi)
library(chipseq)
library(ChIPpeakAnno)
library(GenomicRanges)
library(ggExtra)
library(IRanges)
library(org.Hs.eg.db)
library(reshape2)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

source('code/analysis_functions/custom_plot_theme.R')

###############################################################################################################

select_p53_targets_chip_seq_max_peaks <- function() {
  
  p53_chip <- read.csv('data/chip_seq/hafner_2017_suppl_data_2_p53_chip.csv', head = TRUE)
  p53_transcriptional_targets <- as.character(read.csv('data/rna_seq/hafner2017.csv')$gene)
  p53_chip_p53_targets <- p53_chip[p53_chip$ClosestGene %in% p53_transcriptional_targets,]
  
  p53_chip_p53_targets_max_peak <- data.frame()
  
  for (gene in 1:length(unique(p53_chip_p53_targets$ClosestGene))) {
    
    p53_chip_p53_targets_gene <- 
      p53_chip_p53_targets[p53_chip_p53_targets$ClosestGene == unique(p53_chip_p53_targets$ClosestGene)[gene],]
    p53_chip_p53_targets_max_peak <- 
      rbind(p53_chip_p53_targets_max_peak,
            p53_chip_p53_targets_gene[which.max(p53_chip_p53_targets_gene$NormReads_p53_t0),])
  }
  
  p53_chip_p53_targets_max_peak <- 
    p53_chip_p53_targets_max_peak[order(p53_chip_p53_targets_max_peak$ClosestGene),]
  
  return(p53_chip_p53_targets_max_peak)
}

find_tp53_mdm2_peak_overlaps <- function(cell_line) {
  
  data_directory <- 'data/chip_seq/GSE213300_RAW'
  
  if (cell_line == 'LPS141_untreated') {
    
    tp53_rep1_bed <- import.bed(file.path(data_directory, 'GSM6578472_220505_LPS141_untreated_rep1_P53_S7.bed'))
    tp53_rep2_bed <- import.bed(file.path(data_directory, 'GSM6578476_220505_LPS141_untreated_rep2_P53_S8.bed'))
    
    mdm2_rep1_bed <- import.bed(file.path(data_directory, 'GSM6578471_220505_LPS141_untreated_rep1_MDM2_S1.bed'))
    mdm2_rep2_bed <- import.bed(file.path(data_directory, 'GSM6578475_220505_LPS141_untreated_rep2_MDM2_S2.bed'))
    
    cell_line_name <- 'LPS141'
    
  } else if (cell_line == 'LPS141_nutlin_2h') {
    
    tp53_rep1_bed <- import.bed(file.path(data_directory, 'GSM6578464_220505_LPS141_2h_rep1_P53_S9.bed'))
    tp53_rep2_bed <- import.bed(file.path(data_directory, 'GSM6578468_220505_LPS141_2h_rep2_P53_S10.bed'))
    
    mdm2_rep1_bed <- import.bed(file.path(data_directory, 'GSM6578463_220505_LPS141_2h_rep1_MDM2_S3.bed'))
    mdm2_rep2_bed <- import.bed(file.path(data_directory, 'GSM6578467_220505_LPS141_2h_rep2_MDM2_S4.bed'))
    
    cell_line_name <- 'LPS141'
    
  } else if (cell_line == 'LPS141_nutlin_24h') {
    
    tp53_rep1_bed <- import.bed(file.path(data_directory, 'GSM6578456_220505_LPS141_24h_rep1_P53_S11.bed'))
    tp53_rep2_bed <- import.bed(file.path(data_directory, 'GSM6578460_220505_LPS141_24h_rep2_P53_S12.bed'))
    
    mdm2_rep1_bed <- import.bed(file.path(data_directory, 'GSM6578455_220505_LPS141_24h_rep1_MDM2_S5.bed'))
    mdm2_rep2_bed <- import.bed(file.path(data_directory, 'GSM6578459_220505_LPS141_24h_rep2_MDM2_S6.bed'))
    
    cell_line_name <- 'LPS141'
    
  } else if (cell_line == 'LPS853') {
    
    tp53_rep1_bed <- import.bed(file.path(data_directory, 'GSM6836682_210223_LPS853_DMSO_p53.bed'))
    tp53_rep2_bed <- import.bed(file.path(data_directory, 'GSM6836692_210308_LPS853_DMSO_p53.bed'))
    
    mdm2_rep1_bed <- import.bed(file.path(data_directory, 'GSM6836688_210308_LPS853_DMSO_MDM2.bed'))
    mdm2_rep2_bed <- import.bed(file.path(data_directory, 'GSM6836702_210127_LPS853_200910_DMSO_201208_MDM2.bed'))
    
    cell_line_name <- 'LPS853'
    
  } else if (cell_line == 'HCT116') {
    
    tp53_rep1_bed <- import.bed(file.path(data_directory, 'GSM6836717_220401_HCT116_rep_2_untreated_P53_S18.bed'))
    tp53_rep2_bed <- import.bed(file.path(data_directory, 'GSM6836719_220401_HCT116_rep_3_untreated_P53_S25.bed'))
    
    mdm2_rep1_bed <- import.bed(file.path(data_directory, 'GSM6836713_220317_HCT116_untreated_MDM2_S19.bed'))
    mdm2_rep2_bed <- import.bed(file.path(data_directory, 'GSM6836716_220401_HCT116_rep_2_untreated_MDM2_S17.bed'))
    
    cell_line_name <- 'HCT116'
    
  } else if (cell_line == 'U2OS') {
    
    tp53_rep1_bed <- import.bed(file.path(data_directory, 'GSM6836723_220324_U2OS_untreated_P53_S18.bed'))
    tp53_rep2_bed <- import.bed(file.path(data_directory, 'GSM6836727_220401_U2OS_rep_2_untreated_P53_S22.bed'))
    
    mdm2_rep1_bed <- import.bed(file.path(data_directory, 'GSM6836722_220324_U2OS_untreated_MDM2_S17.bed'))
    mdm2_rep2_bed <- import.bed(file.path(data_directory, 'GSM6836726_220401_U2OS_rep_2_untreated_MDM2_S21.bed'))
    
    cell_line_name <- 'U2OS'
  }
  
  tp53_peaks <- findOverlapsOfPeaks(tp53_rep1_bed, tp53_rep2_bed)$peaklist$`tp53_rep1_bed///tp53_rep2_bed`
  mdm2_peaks <- findOverlapsOfPeaks(mdm2_rep1_bed, mdm2_rep2_bed)$peaklist$`mdm2_rep1_bed///mdm2_rep2_bed`
  tp53_mdm2_overlapping_peaks <- findOverlapsOfPeaks(tp53_peaks, mdm2_peaks)
  
  png(paste0(paste0('figures/', cell_line), '.png'), width = 8, height = 8, units = 'in', res = 120)
  
  makeVennDiagram(
    tp53_mdm2_overlapping_peaks,
    main = cell_line_name,
    main.cex = 2.5,
    main.fontfamily = 'sans',
    main.fontface = 'bold',
    main.pos = c(0.5, 0.95),
    NameOfPeaks = c('TP53', 'MDM2'),
    alpha = 0.75,
    cex = 2.5,
    fontfamily = 'sans',
    fontface = 'bold',
    pos = c(0, 0),
    cat.cex = 2.5,
    cat.fontfamily = 'sans',
    cat.fontface = 'bold',
    cat.pos = c(0, 0),
    fill = c(plot_colors[1], plot_colors[2]))
  
  dev.off()
  
  return(tp53_mdm2_overlapping_peaks)
}

compare_tf_binding_site_to_tss_distances <- function(tp53_mdm2_overlapping_peaks, cell_line_name) {
  
  # Create annotation file from TxDb
  annotation_data <- toGRanges(TxDb.Hsapiens.UCSC.hg38.knownGene, feature = 'gene')
  
  distance_to_tss <- 
    lapply(list(TP53 = tp53_mdm2_overlapping_peaks$peaklist$tp53_peaks,
                MDM2 = tp53_mdm2_overlapping_peaks$peaklist$mdm2_peaks,
                TP53.MDM2 = tp53_mdm2_overlapping_peaks$peaklist$`tp53_peaks///mdm2_peaks`),
           function(.ele) {
             .ele <- annotatePeakInBatch(.ele, AnnotationData = annotation_data, output = 'nearestLocation')
             .ele$shortestDistance})
  distance_to_tss <- melt(distance_to_tss)
  colnames(distance_to_tss) <- c('distance', 'transcription_factor')
  distance_to_tss$transcription_factor <- factor(distance_to_tss$transcription_factor,
                                                 levels = c('TP53', 'MDM2', 'TP53.MDM2'))
  
  p <- ggplot(distance_to_tss, aes(x = transcription_factor, y = distance)) +
    geom_violin(aes(fill = transcription_factor)) +
    geom_boxplot(width = 0.2) +
    scale_y_log10() +
    xlab('Protein') + ylab('Distance to TSS (bp)') + ggtitle(cell_line_name) +
    custom_plot_theme + theme(legend.position = 'none') +
    stat_compare_means(comparisons = list(c('TP53', 'MDM2'), c('MDM2', 'TP53.MDM2'), c('TP53', 'TP53.MDM2')), size = 5)
  
  return(p)
}
