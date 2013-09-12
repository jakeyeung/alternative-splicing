# Jake Yeung
# plot_gene_interactions.R
# 23 June 2013
# Plots gene neighbors from HPRD network.


library(igraph)


# CommandArgs -------------------------------------------------------------

# args <- commandArgs(trailingOnly=TRUE)
# network_fullpath <- args[1]
# genes_of_interest <- args[2]
# plot_output <- args[3]

# Common module btwn two optdis outputs is our genes of interest for 
# Module 3 (gene_exprs_only) and Module 12 (probe_and_gene), respectively
genelist = c('FBXW7', 'UBB', 'SKP2', 'CCNE1')

# FBXW7
# UBB
# SKP2
# CCNE1
# 
# FBXW7
# TFDP1
# RBL1
# CCNE1


# SetDirectories ----------------------------------------------------------

# Set directories for project

# Assumes current directory is source directory
fileDir <- getwd()    # /alternative-splicing/R/
workDir <- dirname(dirname(dirname(fileDir)))    # /alternative-splicing/
inputDir <- file.path(workDir, 'input', 'ppi_networks')
outputDir <- file.path(workDir, 'output', 'plots')



# ReadNetworkFile ---------------------------------------------------------


network <- read.table(file.path(inputDir, 'HPRDNetworkWithComplexes_GeneID.txt'),
                      header=FALSE, sep='\t')
head(network)


# GetGraphObject ----------------------------------------------------------

g <- graph.data.frame(network, directed=FALSE)


# PlotNeighborsofGenes ----------------------------------------------------

for (gene in genelist){
    gene_neighbors <- neighborhood(g, 1, nodes=V(g)[gene])
    sg1 <- induced.subgraph(g, vids=unlist(gene_neighbors))
    print(V(g)[gene_neighbors[[1]]])
    plot(sg1, main=paste(gene, 'Total Neighbors:', length(unlist(gene_neighbors))))
}

gene_names <- c('FBXW7','TFDP1','RBL1','CCNE1')
sg2 <- induced.subgraph(g, vids=gene_names)
plot(sg2)

gene_names <- c('FBXW7', 'UBB', 'SKP2', 'CCNE1')
sg2 <- induced.subgraph(g, vids=gene_names)
plot(sg2)

# FBXW7
# TFDP1
# RBL1
# CCNE1




# gene_names <- c('HRAS', 'VTN', 'PVR', 'TP73', 'INSR', 'IGF2', 'INSRR')
# sg2 <- induced.subgraph(g, vids=gene_names)
# plot(sg2)
# 
# gene_names <- c('GRB7', 'VTN', 'INSR', 'IGF2', 'INSRR', 'HTN1')
# sg2 <- induced.subgraph(g, vids=gene_names)
# plot(sg2)
# 
# gene_names <- c('IL6ST', 'PRKCD', 'PPARA', 'PIK3CA')
# sg2 <- induced.subgraph(g, vids=gene_names)
# plot(sg2)
# 
# gene_names <- c('CREBBP', 'NKX2-1', 'TP73', 'PROX1', 'MYBL1')
# sg2 <- induced.subgraph(g, vids=gene_names)
# plot(sg2)

# CREBBP    21
# NKX2-1	21
# TP73	21
# PROX1	21
# MYBL1	21
# 
# IL6ST    20
# PRKCD	20
# PPARA	20
# PIK3CA	20

# HRAS    3
# VTN	3
# PVR	3
# TP73	3
# INSR	3
# IGF2	3
# INSRR	3

# GRB7    12
# VTN	12
# INSR	12
# IGF2	12
# INSRR	12
