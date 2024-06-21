library(R.matlab)

################################################################################################################

read_rna_seq_data <- function() {
  
  data <- readMat('data/rna_seq/MCF7_Rep1_Rep2.mat')
  genes <- as.character(data.frame(unlist(data$ChuMCFP))[, 1])
  
  mrna_dynamics_p1 <- data.frame(t(data$NGEchMCFP[1:14,]))
  mrna_dynamics_p2 <- data.frame(t(data$NGEchMCFP[35:48,]))
  mrna_dynamics_s1 <- data.frame(t(data$NGEchMCFP[c(1:3, 16:26),]))
  mrna_dynamics_s2 <- data.frame(t(data$NGEchMCFP[c(35:37, 50:60),]))
  
  colnames(mrna_dynamics_p1) <- 
    c('0h', '1h', '2h', '3h', '4h', '5h', '6h', '7h', '8h', '9h', '10h', '11h', '12h', '24h')
  colnames(mrna_dynamics_p2) <- 
    c('0h', '1h', '2h', '3h', '4h', '5h', '6h', '7h', '8h', '9h', '10h', '11h', '12h', '24h')
  colnames(mrna_dynamics_s1) <- 
    c('0h', '1h', '2h', '3h', '4h', '5h', '6h', '7h', '8h', '9h', '10h', '11h', '12h', '24h')
  colnames(mrna_dynamics_s2) <- 
    c('0h', '1h', '2h', '3h', '4h', '5h', '6h', '7h', '8h', '9h', '10h', '11h', '12h', '24h')
  
  mrna_dynamics_p1$gene <- genes
  mrna_dynamics_p2$gene <- genes
  mrna_dynamics_s1$gene <- genes
  mrna_dynamics_s2$gene <- genes
  
  mrna_dynamics_p1[mrna_dynamics_p1 == 0] <- NA
  mrna_dynamics_p2[mrna_dynamics_p2 == 0] <- NA
  mrna_dynamics_s1[mrna_dynamics_s1 == 0] <- NA
  mrna_dynamics_s2[mrna_dynamics_s2 == 0] <- NA
  mrna_dynamics_p1_complete <- mrna_dynamics_p1[complete.cases(mrna_dynamics_p1),]
  mrna_dynamics_p2_complete <- mrna_dynamics_p2[complete.cases(mrna_dynamics_p2),]
  mrna_dynamics_s1_complete <- mrna_dynamics_s1[complete.cases(mrna_dynamics_s1),]
  mrna_dynamics_s2_complete <- mrna_dynamics_s2[complete.cases(mrna_dynamics_s2),]
  
  # Find genes common to all datasets
  common_genes_rna_seq <- 
    Reduce(intersect,
           list(mrna_dynamics_p1_complete$gene,
                mrna_dynamics_p2_complete$gene,
                mrna_dynamics_s1_complete$gene,
                mrna_dynamics_s2_complete$gene))
  
  mrna_dynamics_p1_complete_common_genes <- 
    mrna_dynamics_p1_complete[mrna_dynamics_p1_complete$gene %in% common_genes_rna_seq,]
  mrna_dynamics_p2_complete_common_genes <- 
    mrna_dynamics_p2_complete[mrna_dynamics_p2_complete$gene %in% common_genes_rna_seq,]
  mrna_dynamics_s1_complete_common_genes <- 
    mrna_dynamics_s1_complete[mrna_dynamics_s1_complete$gene %in% common_genes_rna_seq,]
  mrna_dynamics_s2_complete_common_genes <- 
    mrna_dynamics_s2_complete[mrna_dynamics_s2_complete$gene %in% common_genes_rna_seq,]
  
  mrna_dynamics_data <- 
    list(p1 = mrna_dynamics_p1_complete_common_genes, p2 = mrna_dynamics_p2_complete_common_genes,
         s1 = mrna_dynamics_s1_complete_common_genes, s2 = mrna_dynamics_s2_complete_common_genes)
  
  return(mrna_dynamics_data)
}

normalize_rna_seq_data <- function(mrna_dynamics, method = 'percentile') {
  
  rownames(mrna_dynamics) <- mrna_dynamics$gene
  mrna_dynamics <- as.matrix(mrna_dynamics[, 1:13])
  
  if (method == 'none') {
    normalized_mrna_dynamics <- mrna_dynamics
  } else if (method == 'log') {
    normalized_mrna_dynamics <- log2(mrna_dynamics)
  } else if (method == 'fc') {
    normalized_mrna_dynamics <- mrna_dynamics/mrna_dynamics[, 1]
  } else if (method == 'log_fc') {
    normalized_mrna_dynamics <- log2(mrna_dynamics/mrna_dynamics[, 1])
  } else if (method == 'percentile') {
    normalized_mrna_dynamics <- mrna_dynamics/rowMaxs(mrna_dynamics)*100
  } else if (method == 'z_scale') {
    normalized_mrna_dynamics <- (mrna_dynamics - rowMeans(mrna_dynamics))/rowSds(mrna_dynamics) 
  }
  
  return(normalized_mrna_dynamics)
}
