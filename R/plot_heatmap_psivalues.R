# Jake Yeung
# heatmap_psi_values.R
# Heatmaps psi values.

library('RColorBrewer')
library('gplots')
library('Heatplus')

source(file.path(getwd(), 'Rutilities', 'heatmap3.R'))


# Define functions --------------------------------------------------------

set_params_plot_heatmap3 <- function(psi_matrix, plot_title, samp1_name, 
                                     samp2_name){
    #Define custom dist and hclust functions for use with heatmaps
    mydist=function(c) {dist(c,method="euclidian")}
    myclust=function(c) {hclust(c,method="average")}
    main_title <- plot_title
    
    # Init colors
    col_ylgnbu <- colorRampPalette(brewer.pal(n = 9, "YlGnBu"))
    
    # Get subtype colors
    subtype_colors <- as.matrix(data.frame('Subtype'=c(rep('red', ncol(psi_matrix)-1), 
                                                       rep('blue', 1))))
    # pdf(outfile, width=8.5, height=8.5)
    # Plot heatmap3
    heatmap.3(psi_matrix, 
              hclustfun=function(x){hclust(x, method='centroid')}, 
              Colv = FALSE,
              col=col_ylgnbu(256), 
              trace='none',
              margins = c(2,17),
              KeyValueName='PSI Value',
              ColSideColors=subtype_colors,
              NumColSideColors = ncol(subtype_colors),
              labRow='',
              labCol='',
              xlab='Samples',
              density.info = 'density',
              main=main_title)
    legend('topright', legend=c(samp1_name, samp2_name), 
           fill=c('red', 'blue'), 
           border=FALSE, 
           bty='n', 
           y.intersp=0.9, 
           cex=1)
    # dev.off()
}


# Main --------------------------------------------------------------------

fpath.SE <- file.path('G:/jyeung/projects/alternative_splicing/output/',
    'miso_outputs/xenographs_331_331R_hg19_v2/SE.hg19.gff3/',
    'filtered_batch/numinc_2.numexc_2.numsum_5/',
    'tophat_AB331_3_PE_vs_tophat_AB331_R_3_PE.for_heatmap.matrix')
mytitle.SE <- paste0('\nSkipped Exon Events: PSI Values\n')
fpath.MXE <- file.path('G:/jyeung/projects/alternative_splicing/output/',
                   'miso_outputs/xenographs_331_331R_hg19_v2/MXE.hg19.gff3/',
                   'filtered_batch/numinc_2.numexc_2.numsum_5/',
                   'tophat_AB331_3_PE_vs_tophat_AB331_R_3_PE.for_heatmap.matrix')
mytitle.MXE <- paste0('\nMutually Exclusive Events: PSI Values\n')
samp1_name <- 'PCa_331'
samp2_name <- 'NEPC_331R'
lapply(list(c(fpath.SE, mytitle.SE, samp1_name, samp2_name), 
            c(fpath.MXE, mytitle.MXE, samp1_name, samp2_name)),
       function(x){
           # x[1] is filepath, x[2] is title of plot,
           # x[3] is samp1name, x[4] is samp2name
           mymat <- read.table(x[1], header=TRUE, sep='\t')
           genes <- mymat$gsymbol
           psivals <- as.matrix(mymat[, -1])
           plot_title <- paste(x[2], nrow(psivals), 'Events\n', ncol(psivals), 'Samples\n')
           set_params_plot_heatmap3(psivals, plot_title, 
                                    x[3], x[4])
       })

