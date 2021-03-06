# Jake Yeung
# May 27 2013
# parse_dat.R

parse_dat <- function(filename){
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
    
    filename <- filename
    # filename <- 'si_and_gene_expression_data_long.txt'
    unwanted_samples <- c('C42.RNA', 'LN.AI.Luc', 'X890.LN', 'X890.LN', 'X945L.LN')
    
    
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
    
    return(si_gene_data)
}