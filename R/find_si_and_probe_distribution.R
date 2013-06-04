# Jake Yeung
# May 27 2013
# find_si_and_probe_distribution.R
# After running python pipeline, look at the outputs in terms of how the SI
# and gene expression ranges. Perhaps look into ways to modify or
# transform the gene expression to have same dynamic range as SI.

# SetDirectories ----------------------------------------------------------


# Assumes current directory is source directory
fileDir <- getwd()    # /alternative-splicing/R/
sourceDir <- file.path(fileDir, 'Rutilities')
baseDir <- dirname(dirname(dirname(fileDir)))    # /alternative-splicing/
plotDir <- file.path(baseDir, 'output', 'plots')
tableDir <- file.path(baseDir, 'output', 'tables')


# SourceFiles ----------------------------------------------------------

source(file.path(sourceDir, 'parse_and_plot.R'))


# SetConstants ------------------------------------------------------------

filename1 <- 'si_and_gene_expression_data_long_probes.txt'
filename2 <- 'si_and_gene_expression_data_long.txt'

# Output filnames
violin_si_fname <- 'violinplot_si_gene.pdf'
density_si_fname <- 'densityplot_si_gene.pdf'
qqplot_si_fname <- 'qqplot_si_gene.pdf'

violin_probe_gene_fname <- 'violinplot_probe.pdf'
density_probe_gene_fname <- 'densityplot_probe.pdf'
qqplot_probe_fname <- 'qqplot_probe_gene.pdf'

violin_scaled_si_fname <- 'violinplot_scaled_si_gene.pdf'
density_scaled_si_fname <- 'densityplot_scaled_si_gene.pdf'
qqplot_scaled_si_fname <- 'qqplot_scaled_si_gene.pdf'

violin_scaled_probe_fname <- 'violinplot_scaled_probe.pdf'
density_scaled_probe_fname <- 'densityplot_scaled_probe_gene.pdf'
qqplot_scaled_probe_fname <- 'qqplot_scaled_probe_gene.pdf'


# RunFunction -------------------------------------------------------------

# ReadTables
# probe_gene_data <- read.table(file.path(tableDir, filename1), 
#                            header=TRUE, sep='\t')
# si_gene_data <- read.table(file.path(tableDir, filename2), 
#                            header=TRUE, sep='\t')


# Removes unwanted samples. 
probe_gene_data <- parse_dat(filename1)
si_gene_data <- parse_dat(filename2)


# Standardize -------------------------------------------------------------

# Standardize by mean and standard deviation. 
standardized_probes.vector <- tapply(probe_gene_data$probe_or_geneExp, probe_gene_data$gene_symbol, function(x){
    return(scale(x, center=TRUE, scale=TRUE))
})
standardized_probes.matrix <- do.call(rbind, standardized_probes.vector)

standardized_si.vector <- tapply(si_gene_data$probe_or_geneExp, si_gene_data$gene_symbol, function(x){
    return(scale(x, center=TRUE, scale=TRUE))
})
standardized_si.matrix <- do.call(rbind, standardized_si.vector)

# Cbind to probe_gene_data and si_gene_data, respectively.
probe_gene_data <- cbind(probe_gene_data, standardized_probes.matrix)
si_gene_data <- cbind(si_gene_data, standardized_si.matrix)

# Rename columnname
colnames(si_gene_data)[ncol(si_gene_data)] <- 'standardizedExp'
colnames(probe_gene_data)[ncol(probe_gene_data)] <- 'standardizedExp'

# PlotResultsUnstandardized -------------------------------------------------------------


# Compare SI vs Gene Unscaled ---------------------------------------------------

# Violin plot
qplot(factor(probe_or_gene), probe_or_geneExp, data=si_gene_data, geom='violin') + 
    xlab('Gene or SI Value') + 
    ylab('SI Value or Gene Expression') + 
    ggtitle('Violin Plot of SI Value and Gene Expression')
ggsave(file.path(plotDir, violin_si_fname))

qplot(probe_or_geneExp, data=si_gene_data, geom='density') + 
    facet_grid(probe_or_gene ~.) + 
    xlab('Gene or SI Value') + 
    ylab('Density') + 
    ggtitle('Density Plot of SI Value and Gene Expression')
ggsave(file.path(plotDir, density_si_fname))

qplot(sample=probe_or_geneExp, data=si_gene_data, colour=factor(probe_or_gene)) + 
    xlab('Theoretical Quantiles') + 
    ylab('Experimental Quantiles') + 
    ggtitle('QQPlot of SI Values and Gene Expression')
ggsave(file.path(plotDir, qqplot_si_fname))


# Compare Probe vs Gene Unscaled------------------------------------------------------


qplot(factor(probe_or_gene), probe_or_geneExp, data=probe_gene_data, geom='violin') + 
    xlab('Gene or Probe') + 
    ylab('Probe Expression or Gene Expression') + 
    ggtitle('Violin Plot of Probe and Gene Expression')
ggsave(file.path(plotDir, violin_probe_gene_fname))

qplot(probe_or_geneExp, data=probe_gene_data, geom='density') + 
    facet_grid(probe_or_gene ~.) + 
    xlab('Probe Expression or Gene Expression') + 
    ylab('Density') + 
    ggtitle('Density Plot of Probe and Gene Expression')
ggsave(file.path(plotDir, density_probe_gene_fname))

qplot(sample=probe_or_geneExp, data=probe_gene_data, colour=factor(probe_or_gene)) + 
    xlab('Theoretical Quantiles') + 
    ylab('Experimental Quantiles') + 
    ggtitle('QQPlot of Probe and Gene Expression')
ggsave(file.path(plotDir, qqplot_probe_fname))

# PlotResultsStandardized -------------------------------------------------------------


# Compare SI vs Gene Scaled ---------------------------------------------------

qplot(factor(probe_or_gene), standardizedExp, data=si_gene_data, geom='violin') + 
    xlab('Mean-SD Scaled Gene or SI Value') + 
    ylab('SI Value or Gene Expression') + 
    ggtitle('Violin Plot of SI Value and Gene Expression: SD and Mean Scaled')
ggsave(file.path(plotDir, violin_scaled_si_fname))

qplot(standardizedExp, data=si_gene_data, geom='density') + 
    facet_grid(probe_or_gene ~.) + 
    xlab('Mean-SD Scaled Gene or SI Value') + 
    ylab('Density') + 
    ggtitle('Density Plot of SI Value and Gene Expression: SD and Mean Scaled')
ggsave(file.path(plotDir, density_scaled_si_fname))

qplot(sample=standardizedExp, data=si_gene_data, colour=factor(probe_or_gene)) + 
    xlab('Theoretical Quantiles') + 
    ylab('Experimental Quantiles') + 
    ggtitle('QQPlot of SI Values and Gene Expression: SD and Mean Scaled')
ggsave(file.path(plotDir, qqplot_scaled_si_fname))


# Compare Probe vs Gene Scaled ------------------------------------------------------


qplot(factor(probe_or_gene), standardizedExp, data=probe_gene_data, geom='violin') + 
    xlab('Gene or Probe') + 
    ylab('Mean-SD Scaled Probe Expression or Gene Expression') + 
    ggtitle('Violin Plot of Probe and Gene Expression: Mean and SD Scaled')
ggsave(file.path(plotDir, violin_scaled_probe_fname))

qplot(standardizedExp, data=probe_gene_data, geom='density') + 
    facet_grid(probe_or_gene ~.) + 
    xlab('Mean-SD Scaled Gene or Probe') + 
    ylab('Density') + 
    ggtitle('Density Plot of Probe and Gene Expression: Mean and SD Scaled')
ggsave(file.path(plotDir, density_scaled_probe_fname))

qplot(sample=standardizedExp, data=probe_gene_data, colour=factor(probe_or_gene)) + 
    xlab('Theoretical Quantiles') + 
    ylab('Experimental Quantiles') + 
    ggtitle('QQPlot of Mean-SD Scaled SI Values and Gene Expression')
ggsave(file.path(plotDir, qqplot_scaled_probe_fname))

dev.off()

