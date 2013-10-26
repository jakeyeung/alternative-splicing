# Jake Yeung
# heatmap_psi_values.R
# Heatmaps psi values.

library('RColorBrewer')
library('gplots')
library('Heatplus')

source(file.path(getwd(), 'Rutilities', 'heatmap3.R'))


# Define functions --------------------------------------------------------

prepare_t_test_matrix <- function(mydf, not_na_count_threshold=9){
    # Prepare data to be heatplotted.
    # Not NA count threshold: sometimes heatplot won't cluster properly
    # if there are too many NAs for a single gene.
    # We will need to remove those genes for clustering.
    # Remove cols and rows with ALL NAs.
    psi_matrix <- as.matrix(mydf)
    psi_matrix = psi_matrix[rowSums(!is.na(psi_matrix))!=0, 
                            colSums(!is.na(psi_matrix))!=0]
    str(psi_matrix)    # 3 samples have no data.
    
    # Try to remove rows with no variability.
    notnacount <- apply(psi_matrix, 1, function(x){length(which(!is.na(x)))})
    psi_matrix <- psi_matrix[which(notnacount > 9), ]
    return(psi_matrix)
}

set_params_plot_heatmap3 <- function(psi_matrix, plot_title){
    #Define custom dist and hclust functions for use with heatmaps
    mydist=function(c) {dist(c,method="euclidian")}
    myclust=function(c) {hclust(c,method="average")}
    main_title <- plot_title
    
    # Init colors
    col_ylgnbu <- colorRampPalette(brewer.pal(n = 9, "YlGnBu"))
    
    # Get subtype colors
    subtype_colors <- as.matrix(data.frame('Subtype'=c(rep('red', ncol(psi_matrix)-4), 
                                                       rep('blue', 4))))
    
    # Plot heatmap3
    heatmap.3(psi_matrix, 
              hclustfun=function(x){hclust(x, method='centroid')}, 
              Colv = FALSE,
              col=col_ylgnbu(256), 
              trace='none',
              margins = c(4,17),
              KeyValueName='PSI Value',
              ColSideColors=subtype_colors,
              NumColSideColors = ncol(subtype_colors),
              labRow='',
              labCol='',
              xlab='Patients',
              density.info = 'histogram',
              main=main_title)
    legend('topright', legend=c('Adenocarcinoma', 'Neuroendocrine'), 
           fill=c('red', 'blue'), 
           border=FALSE, 
           bty='n', 
           y.intersp=0.9, 
           cex=1)
}


# Main --------------------------------------------------------------------


matrix_file.SE <- paste0('G:/jyeung/projects/',
                         'alternative_splicing/output/miso_outputs/',
                         'mark_rubin_hg19_v2_rl_insertdist/SE.hg19.gff3/',
                         't_test_results/t_test_mincounts_10.matrix')

matrix_file.MXE <- paste0('G:/jyeung/projects/',
                          'alternative_splicing/output/miso_outputs/',
                          'mark_rubin_hg19_v2_rl_insertdist/MXE.hg19.gff3/',
                          't_test_results/t_test_mincounts_10.matrix')

se_plot_title <- 'PSI-Values Across Prostate Cancer Patients:\nSkipped Exon Events'
mxe_plot_title <- 'PSI-Values Across Prostate Cancer Patients:\nMutually Exclusive Events'

lapply(list(c(matrix_file.SE, se_plot_title), c(matrix_file.MXE, mxe_plot_title)), 
       function(x){
           # x[1] is filepath, x[2] is title of plot
        mymat <- as.matrix(read.table(x[1], row.names=1, header=TRUE, sep='\t'))
        mymat <- prepare_t_test_matrix(mymat)
        set_params_plot_heatmap3(mymat, x[2])
})
# psi_matrix.SE <- as.matrix(read.table(matrix_file, row.names=1, header=TRUE, sep='\t'))
# psi_matrix.SE <- prepare_t_test_matrix(psi_matrix)
# set_params_plot_heatmap3(psi_matrix)
# 
# psi_matrix.MXE <- as.matrix(read.table(matrix_file, row.names=1, header=TRUE, sep='\t'))
# psi_matrix.MXE <- prepare_t_test_matrix(psi_matrix)
# set_params_plot_heatmap3(psi_matrix)



