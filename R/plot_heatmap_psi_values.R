# Jake Yeung
# plot_heatmap_psi_values.R
# Heatmaps psi values.
# June 23 2014
# usage: Rscript plot_heatmap_psi_values.R matrix.txt outputfile.eps

library('RColorBrewer')
library('gplots')
library("extrafont")    # for Arial font

source(file.path(getwd(), 'Rutilities', 'heatmap3.R'))


# Define functions --------------------------------------------------------

set_params_plot_heatmap3 <- function(psi_matrix, plot_title, samp1_name, 
                                     samp2_name){
    #Define custom dist and hclust functions for use with heatmaps
    mydist=function(c) {dist(c,method="euclidian")}
    myclust=function(c) {hclust(c,method="average")}
    main_title <- plot_title
    
    # Init colors
    col_ylgnbu <- colorRampPalette(rev(brewer.pal(n = 9, "YlGnBu")))
    
    # Get subtype colors
    subtype_colors <- as.matrix(data.frame('Subtype'=c(rep('red', ncol(psi_matrix)-1), 
                                                       rep('blue', 1))))
    #subtype_colors <- as.matrix(data.frame('Subtype'=c(rep('blue', ncol(psi_matrix)-1), 
    #rep('red', 1))))
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
           y.intersp=0.7, 
           cex=1)
    # dev.off()
}


# Main -----------------------------------------------------------------    ---

args <- commandArgs(trailingOnly=TRUE)
myfilepath <- args[1]
se_outpath <- args[2]

fpath.SE <- file.path(myfilepath)
mytitle.SE <- paste0('\nCassette Exons\n')
samp1_name <- 'Sample 1'
samp2_name <- 'Sample 2'

loadfonts(device='postscript')    # make things arial

lapply(list(c(fpath.SE, mytitle.SE, samp1_name, samp2_name, se_outpath)),
       function(x){
           # x[1] is filepath, x[2] is title of plot,
           # x[3] is samp1name, x[4] is samp2name
           postscript(x[5], height=9.19, width=6.83,
                      family='Arial', paper='special',
                      onefile=FALSE,
                      horizontal=FALSE)
           mymat <- read.table(x[1], header=TRUE, sep='\t')
           genes <- mymat$gsymbol
           psivals <- as.matrix(mymat[, -1])
           plot_title <- paste('\n', x[2], nrow(psivals), 'Events\n', ncol(psivals), 'Samples\n')
           set_params_plot_heatmap3(psivals, plot_title, 
                                    x[3], x[4])
           dev.off()
           print(paste(c('Figure output:', x[5])))
       })


