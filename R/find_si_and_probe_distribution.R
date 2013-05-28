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


# RunFunction -------------------------------------------------------------

probe_gene_data <- parse_dat(filename1)
si_gene_data <- parse_dat(filename2)


# PlotResults -------------------------------------------------------------


# Compare Probe vs Gene ---------------------------------------------------

qplot(factor(probe_or_gene), probe_or_geneExp, data=si_gene_data, geom='violin')

qplot(probe_or_geneExp, data=si_gene_data, geom='density') + 
    facet_grid(probe_or_gene ~.)

qplot(sample=probe_or_geneExp, data=si_gene_data, colour=factor(probe_or_gene))


# Compare SI vs Gene ------------------------------------------------------


qplot(factor(probe_or_gene), probe_or_geneExp, data=probe_gene_data, geom='violin')

qplot(probe_or_geneExp, data=probe_gene_data, geom='density') + 
    facet_grid(probe_or_gene ~.)

qplot(sample=probe_or_geneExp, data=probe_gene_data, colour=factor(probe_or_gene))
