# Jake Yeung
# heatmap_psi_values.R
# Heatmaps psi values.

library('RColorBrewer')
library('gplots')
#library('Heatplus')
library("extrafont")    # for Arial font

source(file.path(getwd(), 'Rutilities', 'heatmap3.R'))


# Define functions --------------------------------------------------------

# NA values are replaced
# recommended by:
# https://stat.ethz.ch/pipermail/bioconductor/2003-October/002874.html
na.dist <- function(x,...) {
    t.dist <- dist(x,...)
    t.dist <- as.matrix(t.dist)
    t.limit <- 1.1*max(t.dist,na.rm=T)
    t.dist[is.na(t.dist)] <- t.limit
    t.dist <- as.dist(t.dist)
    return(t.dist)
}

prepare_t_test_matrix <- function(mydf, not_na_count_threshold=9){
    # Prepare data to be heatplotted.
    # Not NA count threshold: sometimes heatplot won't cluster properly
    # if there are too many NAs for a single gene.
    # We will need to remove those genes for clustering.
    # Remove cols and rows with ALL NAs.
    psi_matrix <- as.matrix(mydf)
    psi_matrix = psi_matrix[rowSums(!is.na(psi_matrix))!=0, 
                            colSums(!is.na(psi_matrix))!=0]
    # str(psi_matrix)    # 3 samples have no data.
    
    # Try to remove rows with no variability.
    notnacount <- apply(psi_matrix, 1, function(x){length(which(!is.na(x)))})
    psi_matrix <- psi_matrix[which(notnacount > not_na_count_threshold), ]
    return(psi_matrix)
}

set_params_plot_heatmap3 <- function(psi_matrix, plot_title, n_neuro_samps=n_neuro_samps){
    #Define custom dist and hclust functions for use with heatmaps
    mydist=function(c) {dist(c,method="euclidian")}
    myclust=function(c) {hclust(c,method="average")}
    main_title <- plot_title
    
    # Init colors
    col_ylgnbu <- colorRampPalette(rev(brewer.pal(n = 9, "YlGnBu")))
    
    # Get subtype colors
    subtype_colors <- as.matrix(data.frame('Subtype'=c(rep('red', ncol(psi_matrix)-n_neuro_samps), 
                                                       rep('blue', n_neuro_samps))))
    
    # Plot heatmap3
    heatmap.3(psi_matrix,
              #revC = TRUE,
              hclustfun=function(x){hclust(na.dist(x), method='centroid')}, 
              Colv = FALSE,
              col=col_ylgnbu(256), 
              trace='none',
              margins = c(4,17),    # margins for no sample names
              #margins = c(12,17),    # if you want to see sample names
              KeyValueName='PSI Value',
              ColSideColors=subtype_colors,
              NumColSideColors = ncol(subtype_colors),
              labRow='',
              labCol='',    # comment out to see sample names
              xlab='Patient Tumour Samples',
              density.info = 'density',
              main=main_title)
    legend('topright', legend=c('Adenocarcinoma', 'Neuroendocrine'), 
           fill=c('red', 'blue'), 
           border=FALSE, 
           bty='n', 
           y.intersp=0.9, 
           cex=1)
}


# Main --------------------------------------------------------------------
matrix_file.SE <- 'G:/jyeung/projects/alternative_splicing/output/miso_outputs/rubin_takeda_pooled/SE.hg19.gff3/t_test_results/filter030/t_test_mincounts_10.pval_adj.exprs_filtered.genenames.filtered_0.30.for_heatmap.matrix'
    
#matrix_file.MXE <- paste0('G:/jyeung/projects/',
#                          'alternative_splicing/output/miso_outputs/',
#                          'rubin_takeda_pooled/MXE.hg19.gff3/',
#                          't_test_results/t_test_mincounts_10.matrix_3_removed.txt')

se_plot_title <- '\n\nCassette Exons\n'
#mxe_plot_title <- '\n\nMutually Exclusive Exons\n'

se_outpath <- 'G:/jyeung/projects/alternative_splicing/figs/cassette_exons_pooled.rerun3.dendrogram.eps'
#mxe_outpath <- 'G:/jyeung/projects/alternative_splicing/figs/mutually_exclusive_pooled_3_removed_rev.eps' 

loadfonts(device='postscript')    # make things arial

lapply(list(c(matrix_file.SE, se_plot_title, se_outpath)), 
       function(x){
           # x[1] is filepath, x[2] is title of plot
        postscript(x[3], height=9.19, width=6.83,
                   family='Arial', paper='special',
                   onefile=FALSE,
                   horizontal=FALSE)
        #op <- par(mar = c(5, 4, 0.05, 0.05) + 0.1)
        mymat <- as.matrix(read.table(x[1], row.names=1, header=TRUE, sep='\t'))
        #genes <- mymat$gsymbol
        #psivals <- as.matrix(mymat[, -1])
        plot_title <- paste(x[2], nrow(mymat), 'Events\n', ncol(mymat), 'Samples\n')
        set_params_plot_heatmap3(mymat, plot_title, n_neuro_samps=9)
        dev.off()
        print(paste(c('Figure output:', x[3])))
})

