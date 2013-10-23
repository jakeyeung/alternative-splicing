# Jake Yeung
# heatmap_psi_values.R
# Heatmaps psi values.

library('RColorBrewer')
library('gplots')


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
ind <- apply(psi_matrix, 1, function(x){var(x, na.rm=TRUE)})
# psi_matrix <- psi_matrix[!ind, ]

# Init colors
col_ylgnbu <- colorRampPalette(brewer.pal(n = 9, "YlGnBu"))
jGraysFun <- colorRampPalette(brewer.pal(n = 9, "Greys"))
heatmap(psi_matrix, Rowv=NA, Colv=NA, col_ylgnbu)

