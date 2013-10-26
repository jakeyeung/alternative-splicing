# Jake Yeung
# heatmap_psi_values.R
# Heatmaps psi values.

library('RColorBrewer')
library('gplots')
library('Heatplus')

source(file.path(getwd(), 'Rutilities', 'heatmap3.R'))


matrix_file <- paste0('G:/jyeung/projects/',
    'alternative_splicing/output/miso_outputs/',
    'mark_rubin_hg19_v2_rl_insertdist/SE.hg19.gff3/',
    't_test_results/t_test_mincounts_10.matrix')

psi_matrix <- as.matrix(read.table(matrix_file, row.names=1, header=TRUE, sep='\t'))
str(psi_matrix)

# Remove cols and rows with ALL NAs.
psi_matrix = psi_matrix[rowSums(!is.na(psi_matrix))!=0, 
                        colSums(!is.na(psi_matrix))!=0]
str(psi_matrix)    # 3 samples have no data.

# Try to remove rows with no variability.
notnacount <- apply(psi_matrix, 1, function(x){length(which(!is.na(x)))})
# nacount_col <- apply(psi_matrix, 2, function(x){length(which(is.na(x)))})
# ind <- apply(psi_matrix, 1, function(x){var(x, na.rm=TRUE)})
# ind <- apply(psi_matrix, 1, var)
psi_matrix <- psi_matrix[which(notnacount > 9), ]

#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}
main_title <- 'PSI-Values Across Prostate Cancer Patients'

# Init colors
col_ylgnbu <- colorRampPalette(brewer.pal(n = 9, "YlGnBu"))
jGraysFun <- colorRampPalette(brewer.pal(n = 9, "Greys"))

# # 
# reg4 <- regHeatmap(psi_matrix, legend=3, dendrogram=list(Col=list(status='hide')), col=col_ylgnbu, scale='none')
# plot(reg4)


# # Heatmap3
# heatmap.3(psi_matrix, 
#           hclustfun=myclust, 
#           distfun=mydist, 
#           na.rm = FALSE, 
#           scale="none", 
#           dendrogram="both", 
#           margins=c(4,10), 
#           Rowv=TRUE, 
#           Colv=TRUE, 
#           ColSideColors=clab,
#           RowSideColors=rlab, 
#           symbreaks=FALSE, 
#           key=TRUE, symkey=FALSE, 
#           density.info="none", 
#           trace="none", 
#           main=main_title, 
#           labCol=FALSE, 
#           labRow=FALSE, 
#           cexRow=1, 
#           col=col_ylgnbu(75), 
#           NumColSideColors=7, 
#           KeyValueName="PSI-Values")

# Init colors
col_ylgnbu <- colorRampPalette(brewer.pal(n = 9, "YlGnBu"))
jGraysFun <- colorRampPalette(brewer.pal(n = 9, "Greys"))

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
          main='PSI-Values Across Prostate Cancer Patients')
legend('topright', legend=c('Adenocarcinoma', 'Neuroendocrine'), 
       fill=c('red', 'blue'), 
       border=FALSE, 
       bty='n', 
       y.intersp=0.9, 
       cex=1)

