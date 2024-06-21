library(circlize)
library(clusterProfiler)
library(cogena)
library(ComplexHeatmap)
library(DOSE)
library(enrichplot)
library(dplyr)
library(GSEABase)
library(GSVA)
library(msigdbr)
library(pathview)
library(ReactomePA)
library(singscore)

################################################################################################################

compute_geneset_variation_analysis_scores <- function(molecular_data, genesets_type, genesets_name, 
                                                      normalize = FALSE, timecourse = FALSE) {
  
  # Retrieve the normalized data from the `DESeqDataSet`
  vst_df <- assay(molecular_data) %>%
    as.data.frame() %>%
    tibble::rownames_to_column('ensembl_id') # Make Gene IDs into their own column
  
  if (genesets_type == 'msigdb') {
    
    if (genesets_name == 'hallmark') {
      
      hallmark_gene_sets <- msigdbr::msigdbr(
        species = 'Homo sapiens',
        category = 'H'
      )
      
      geneset_list <- split(
        hallmark_gene_sets$entrez_gene, # The genes to split into pathways
        hallmark_gene_sets$gs_name # The pathways made as the higher levels of the list
      )
    }
    
    gene_identifier <- 'ENTREZID'
    
  } else if (genesets_type == 'custom') {
    
    geneset_list <- gmt2list(paste0(paste0('data/signatures/', genesets_name), '.gmt'))
    
    gene_identifier <- 'SYMBOL'
  }
  
  vst_df$ensembl_id <- substr(vst_df$ensembl_id, 1, 15)
  
  # Create a mapped data frame to join to the gene expression values
  mapped_df <- data.frame(
    'entrez_id' = mapIds(
      org.Hs.eg.db,
      keys = vst_df$ensembl_id,
      # Replace with the type of gene identifiers in your data
      keytype = 'ENSEMBL',
      # Replace with the type of gene identifiers you would like to map to
      column = gene_identifier,
      # This will keep only the first mapped value for each Ensembl ID
      multiVals = 'first'
    )
  ) %>%
    # If an Ensembl gene identifier does not map to a Entrez gene identifier, drop that from the data frame
    dplyr::filter(!is.na(entrez_id)) %>%
    # Make an `Ensembl` column to store the row names
    tibble::rownames_to_column('Ensembl') %>%
    # Join the rest of the expression data
    dplyr::inner_join(vst_df, by = c('Ensembl' = 'ensembl_id'))
  
  # Determine the gene means
  gene_means <- rowMeans(mapped_df %>% dplyr::select(-Ensembl, -entrez_id))
  
  # Add this as a column in `mapped_df`.
  mapped_df <- mapped_df %>%
    # Add gene_means as a column called gene_means
    dplyr::mutate(gene_means) %>%
    # Reorder the columns so `gene_means` column is upfront
    dplyr::select(Ensembl, entrez_id, gene_means, dplyr::everything())
  
  filtered_mapped_df <- mapped_df %>%
    # Sort so that the highest mean expression values are at the top
    dplyr::arrange(dplyr::desc(gene_means)) %>%
    # Filter out the duplicated rows using `dplyr::distinct()`
    dplyr::distinct(entrez_id, .keep_all = TRUE)
  
  filtered_mapped_matrix <- filtered_mapped_df %>%
    # GSVA cannot use the Ensembl IDs so drop this column as well as the means
    dplyr::select(-Ensembl, -gene_means) %>%
    # Store gene identifiers as row names
    tibble::column_to_rownames('entrez_id') %>%
    # Convert object into a matrix
    as.matrix()
  
  gsva_results <- gsva(
    filtered_mapped_matrix,
    geneset_list,
    method = 'gsva',
    # Appropriate for VST transformed data
    kcdf = 'Gaussian',
    # Minimum gene set size
    min.sz = 3,
    # Maximum gene set size
    max.sz = 500,
    # Compute Gaussian-distributed scores
    mx.diff = TRUE,
    verbose = FALSE
  )
  
  return(gsva_results)
}

plot_geneset_variation_analysis_heatmap <- function(molecular_data, genesets_type, genesets_name, 
                                                    normalize = FALSE, timecourse = FALSE) {
  
  gsva_results <- compute_geneset_variation_analysis_scores(molecular_data, genesets_type, genesets_name)
  
  conditions <- as.data.frame(colData(molecular_data)[, c('treatment', 'time', 'replicate')])
  
  if (normalize == TRUE) {
    gsva_results <- scale(gsva_results, center = TRUE, scale = TRUE) 
  }
  
  if (timecourse == TRUE) {
    column_split <- c(rep('IR 1', ncol(gsva_results)/4),
                      rep('IR 2', ncol(gsva_results)/4),
                      rep('IR + nutlin-3a 1', ncol(gsva_results)/4),
                      rep('IR + nutlin-3a 2', ncol(gsva_results)/4))
  } else {
    column_split <- NULL
  }
  
  column_annotation <- 
    columnAnnotation(
      Treatment = conditions$treatment,
      Time = conditions$time,
      Replicate = conditions$replicate,
      col = list(
        Treatment = c('IR' = plot_colors[1], 'IR + nutlin' = plot_colors[2]),
        Time = colorRamp2(c(3, 12), c('white', 'pink')),
        Replicate = c('1' = 'white', '2' = 'gray')),
      border = TRUE,
      annotation_name_gp = gpar(fontsize = 16),
      annotation_legend_param = list(
        Time = list(title = 'Time (h)'),
        title_gp = gpar(fontsize = 16),
        labels_gp = gpar(fontsize = 16)))
  
  colnames(gsva_results) <- as.character(rep(seq(3, 12), 4))
  
  gsva_heatmap <- 
    Heatmap(
      gsva_results,
      name = 'GSVA',
      col = colorRampPalette(rev(c('#D73027', '#FC8D59', '#FEE090', '#FFFFBF', '#E0F3F8', '#91BFDB', '#4575B4')))(100),
      rect_gp = gpar(col = 'black', lwd = 1),
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      show_column_names = TRUE,
      column_names_rot = 90,
      column_title = 'Time (h)',
      column_title_side = 'bottom',
      top_annotation = column_annotation,
      column_split = factor(c(rep('IR 1', 10), rep('IR 2', 10), rep('IR + nutlin 1', 10), rep('IR + nutlin 2', 10)),
                            levels = c('IR 1', 'IR 2', 'IR + nutlin 1', 'IR + nutlin 2')),
      column_gap = unit(0.025, 'npc'),
      border = TRUE,
      row_names_max_width = unit(9, 'cm'),
      column_title_gp = gpar(fontsize = 16),
      row_title_gp = gpar(fontsize = 16),
      column_names_gp = gpar(fontsize = 16),
      row_names_gp = gpar(fontsize = 16),
      heatmap_legend_param = list(
        title_gp = gpar(fontsize = 16),
        labels_gp = gpar(fontsize = 16)))
  
  return(gsva_heatmap)
}

compute_geneset_variation_analysis_scores_metabric <- function(genesets_name = 'ir_nutlin_vs_ir_signature') {
  
  metabric_expression <- read_tsv('data/radiosensitivity/brca_metabric/data_mrna_illumina_microarray.txt', comment = '#')
  
  metabric_expression <- metabric_expression %>%
    dplyr::select(-'Entrez_Gene_Id')
  
  colnames(metabric_expression)[1] <- 'Gene'
  
  mapped_df <- metabric_expression
  
  gene_means <- rowMeans(mapped_df %>% dplyr::select(-Gene))
  
  # Let's add this as a column in our `mapped_df`.
  mapped_df <- mapped_df %>%
    # Add gene_means as a column called gene_means
    dplyr::mutate(gene_means) %>%
    # Reorder the columns so `gene_means` column is upfront
    dplyr::select(Gene, gene_means, dplyr::everything())
  
  filtered_mapped_df <- mapped_df %>%
    # Sort so that the highest mean expression values are at the top
    dplyr::arrange(dplyr::desc(gene_means)) %>%
    # Filter out the duplicated rows using `dplyr::distinct()`
    dplyr::distinct(Gene, .keep_all = TRUE)
  
  filtered_mapped_matrix <- filtered_mapped_df %>%
    # GSVA can't the Ensembl IDs so we should drop this column as well as the means
    dplyr::select(-gene_means) %>%
    # We need to store our gene identifiers as row names
    tibble::column_to_rownames('Gene') %>%
    # Now we can convert our object into a matrix
    as.matrix()
  
  # Gene set variation analysis
  geneset_list <- gmt2list(paste0(paste0('data/signatures/', genesets_name), '.gmt'))
  
  gsva_results <- gsva(
    filtered_mapped_matrix,
    geneset_list,
    method = 'gsva',
    kcdf = 'Gaussian',
    # Minimum gene set size
    min.sz = 5,
    # Maximum gene set size
    max.sz = 2000,
    # Compute Gaussian-distributed scores
    mx.diff = TRUE,
    verbose = FALSE
  )
  
  gsva_results <- data.frame(t(gsva_results))
  gsva_results$Sample_ID <- row.names(gsva_results)
  
  return(gsva_results)
}