library('AnnotationDbi')
library('DEGreport')
library('DESeq2')
library('dplyr')
library('EnhancedVolcano')
library('genefilter')
library('ggbeeswarm')
library('ggplot2')
library('glmpca')
library('Gviz')
library('IHW')
library('magrittr')
library('org.Hs.eg.db')
library('pheatmap')
library('PoiClaClu')
library('RColorBrewer')
library('ReportingTools')
library(reshape2)
library('stringr')
library(tibble)
library('topGO')
library('tximeta')
library('vsn')

source('code/analysis_functions/custom_plot_theme.R')

################################################################################################################

make_linked_transcriptome <- function() {
  
  # Point to the Salmon index to create a linkedTxome as the index will not match a known txome
  index_dir <- file.path(dir, 'gencode.v33_salmon_1.1.0_index')
  
  # Point to the source FASTA and GTF
  fasta_FTP <- file.path(dir, 'gencode.v33.pc_transcripts.fa.gz')
  gtf_path <- file.path(dir, 'gencode.v33.annotation.gtf')
  
  # now create a linkedTxome, linking the Salmon index to its FASTA and GTF sources
  makeLinkedTxome(indexDir = index_dir, source = 'GENCODE', organism = 'Homo sapiens',
                  release = '33', genome = 'GRCh38', fasta = fasta_FTP, gtf = gtf_path, write = FALSE)
}

create_differential_expression_dataset <- function(coldata, comparison) {
  
  coldata$replicate <- as.factor(coldata$replicate)
  
  summarized_experiment <- tximeta(coldata, importer = read.delim)
  # Summarize the transcript-level quantifications to the gene level
  gene_level_summarized_experiment <- summarizeToGene(summarized_experiment)
  
  gene_level_summarized_experiment$treatment <- factor(gene_level_summarized_experiment$treatment,
                                                       levels = c('IR', 'IR + nutlin'))
  levels(gene_level_summarized_experiment$replicate) <- c('1', '2')
  
  if (comparison == 'ir_nutlin_vs_ir_timecourse') {
    
    differential_expression_dataset <- 
      DESeqDataSet(gene_level_summarized_experiment, design = ~ time + replicate + treatment + treatment:time)
    
  } else if (comparison == 'ir_nutlin_vs_ir_single_timepoint') {
    
    gene_level_summarized_experiment$time <- as.factor(gene_level_summarized_experiment$time)
    
    differential_expression_dataset <- DESeqDataSet(gene_level_summarized_experiment, design = ~ treatment)
    
  } else if (comparison == 'treated_vs_control_timecourse') {
    
    differential_expression_dataset <- 
      DESeqDataSet(gene_level_summarized_experiment, design = ~ time + replicate)
  
  } else if (comparison == 'treated_vs_control_single_timepoint') {
    
    gene_level_summarized_experiment$time <- as.factor(gene_level_summarized_experiment$time)
    
    differential_expression_dataset <- 
      DESeqDataSet(gene_level_summarized_experiment, design = ~ replicate + time)
  }
  
  return(differential_expression_dataset)
}

calculate_differential_expression_timecourse <- 
  function(differential_expression_timecourse_dataset, comparison) {
  # Differential expression analysis for time course experiments
  
  if (comparison == 'ir_nutlin_vs_ir_timecourse') {
    
    differential_expression_timecourse_dataset <-
      DESeq(differential_expression_timecourse_dataset, test = 'LRT', reduced = ~ time + replicate + treatment)
    
  } else if (comparison == 'treated_vs_control_timecourse') {
    
    differential_expression_timecourse_dataset <-
      DESeq(differential_expression_timecourse_dataset, test = 'LRT', reduced = ~ replicate)
  }
  
  resTC <- results(differential_expression_timecourse_dataset)
  resTC$symbol <- mcols(differential_expression_timecourse_dataset)$symbol
  head(resTC[order(resTC$padj),], 10)
  
  print(resultsNames(differential_expression_timecourse_dataset))

  return(resTC)
}

plot_mrna_timecourses <- function(dataset, results) {
  
  mrna_timecourse <- 
    plotCounts(dataset, which.min(results$padj),
               intgroup = c('time', 'treatment'), returnData = TRUE)

  mrna_timecourse$time <- as.numeric(as.character(mrna_timecourse$time))
  p <- ggplot(mrna_timecourse, aes(x = time, y = count, color = treatment, group = treatment)) + 
    geom_point() + stat_summary(fun = mean, geom = 'line') +
    scale_y_log10() +
    custom_plot_theme
  print(p)
}

annotate_results <- function(results) {
  
  columns(org.Hs.eg.db)
  ens.str <- substr(rownames(results), 1, 15)
  results$symbol <- mapIds(org.Hs.eg.db,
                           keys = ens.str,
                           column = 'SYMBOL',
                           keytype = 'ENSEMBL',
                           multiVals = 'first')
  results$entrez <- mapIds(org.Hs.eg.db,
                           keys = ens.str,
                           column = 'ENTREZID',
                           keytype = 'ENSEMBL',
                           multiVals = 'first')
  results_ordered <- results[order(results$pvalue),]
  print(head(results_ordered))
  
  return(results_ordered)
}

export_results <- function(results) {
  
  results_ordered <- annotate_results(results)
  
  results_ordered_df <- as.data.frame(results_ordered)
  write.csv(results_ordered_df, file = 'outputs/rna_seq/differential_expression_results.csv')
  
  htmlRep <- HTMLReport(shortName = 'report', title = 'My report',
                        reportDirectory = 'outputs/rna_seq/differential_expression_report')
  publish(results_ordered_df, htmlRep)
  url <- finish(htmlRep)
  browseURL(url)
}

plot_fold_change_in_genomic_space <- function(differential_expression_dataset) {
  
  resGR <- lfcShrink(differential_expression_dataset,
                     coef = 'treatment_IR_vs_IR + nutlin', type = 'apeglm', format = 'GRanges')
  resGR
  ens.str <- substr(names(resGR), 1, 15)
  resGR$symbol <- mapIds(org.Hs.eg.db, ens.str, 'SYMBOL', 'ENSEMBL')
  window <- resGR[topGene] + 1e6
  strand(window) <- '*'
  resGRsub <- resGR[resGR %over% window]
  naOrDup <- is.na(resGRsub$symbol) | duplicated(resGRsub$symbol)
  resGRsub$group <- ifelse(naOrDup, names(resGRsub), resGRsub$symbol)
  status <- factor(ifelse(resGRsub$padj < 0.05 & !is.na(resGRsub$padj), 'sig', 'notsig'))
  options(ucscChromosomeNames = FALSE)
  g <- GenomeAxisTrack()
  a <- AnnotationTrack(resGRsub, name = 'gene ranges', feature = status)
  d <- DataTrack(
    resGRsub, data = 'log2FoldChange', baseline = 0, type = 'h', name = 'log2 fold change', strand = '+')
  plotTracks(list(g, d, a), groupAnnotation = 'group', notsig = 'grey', sig = 'hotpink')  
}

perform_differential_expression_analysis <- function(coldata, comparison) {
  
  differential_expression_dataset <- create_differential_expression_dataset(coldata, comparison)
  
  if (comparison == 'ir_nutlin_vs_ir_timecourse' |
      comparison == 'treated_vs_control_timecourse') {
    
    differential_expression_results <- 
      calculate_differential_expression_timecourse(differential_expression_dataset, comparison)
    
  } else if (comparison == 'ir_nutlin_vs_ir_single_timepoint' |
             comparison == 'treated_vs_control_single_timepoint') {
    
    differential_expression_results <- 
      calculate_differential_expression(differential_expression_dataset, comparison)
  }
  
  differential_expression_results <- annotate_results(differential_expression_results)
  export_results(differential_expression_results)
  plot_mrna_timecourses(differential_expression_dataset, differential_expression_results)
  
  return(differential_expression_results)
}

plot_volcano_plot_rna_seq_timecourse <- function(differential_expression_results) {
  
  differential_expression_plot_data <- data.frame(differential_expression_results)
  
  p <- 
    ggplot(differential_expression_plot_data, aes(log2FoldChange, -log10(padj), label = symbol)) +
    geom_point(color = 'gray') + 
    geom_point(data = subset(differential_expression_plot_data, log2FoldChange < -0.1 & padj < 0.05),
               color = plot_colors[1]) +
    geom_point(data = subset(differential_expression_plot_data, log2FoldChange > 0.1 & padj < 0.05),
               color = plot_colors[2]) +
    geom_text_repel(data = subset(differential_expression_plot_data,
                                  (log2FoldChange > 0.1 | log2FoldChange < -0.1) & padj < 0.05),
                    max.overlaps = 10) +
    xlim(c(-2, 2)) +
    xlab('Effect size') + ylab(expression(-log[10](p[adj])))
  
  return(p)
}
