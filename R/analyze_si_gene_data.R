# Jake Yeung
# May 23 2013
# analyze_si_gene_data.R
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


# LoadLibraries -----------------------------------------------------------


library(ggplot2)



# SourceFiles ----------------------------------------------------------



# SetConstants ------------------------------------------------------------

filename <- 'si_and_gene_expression_data_long_probes.txt'
# filename <- 'si_and_gene_expression_data_long.txt'
unwanted_samples <- c('C42.RNA', 'LN.AI.Luc', 'X1005', 
                     'X890L', 'X945', 'X961')


# ReadTable ---------------------------------------------------------------


si_gene_data <- read.table(file.path(tableDir, filename), 
                           header=TRUE, sep='\t')
str(si_gene_data)



# ModifyTable -------------------------------------------------------------

# BEGIN: Remove NA gene_symbols
# 1. Get row indices containing NA
na_rowindices <- which(is.na(si_gene_data$gene_symbol)==TRUE)
print(paste(length(na_rowindices), 
            'rows containing NA gene symbols'))    # 5551 rows
# 2. Remove rows according to row indices
si_gene_data <- si_gene_data[-na_rowindices, ]
str(si_gene_data)
# END: Remove NA gene_symbols

# BEGIN: Remove unwanted sample IDs
unwanted_sample_rowindices <- which(si_gene_data$sampleID %in% unwanted_samples)
print(paste(length(unwanted_sample_rowindices), 
      'rows containing unwanted samples'))    # 69402 rows
si_gene_data <- si_gene_data[-unwanted_sample_rowindices, ]
str(si_gene_data)


# PlotResults -------------------------------------------------------------

qplot(factor(probe_or_gene), probe_or_geneExp, data=si_gene_data, geom='violin')

qplot(probe_or_geneExp, data=si_gene_data, geom='density') + 
    facet_grid(probe_or_gene ~.)

qplot(sample=probe_or_geneExp, data=si_gene_data, colour=factor(probe_or_gene))



