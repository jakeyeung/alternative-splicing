'''
Created on 2013-06-06

@author: jyeung
'''


import sys
import os
import csv
from utilities import set_directories, plots


# Set constants
inputdir = 'input'
outputdir = 'output'

gene_symbol_colname = 'gene_symbol'    # Colname in both optdis input and output
probe_or_gene_colname = 'probe_or_gene'    # Colname in OptDis Input
probe_or_gene_colname_optdisresults = 'gene'    # Colname in OptDis Output

barplot_categories = ['Non-AS Genes', 'AS-Genes']
barplot_events = ['Input', 'PPI', 'OptDis Output']
barplot_barwidths = 0.35
color_vector = ['red', 'blue']
xlabel = ''
ylabel_count = 'Number of Genes'
title_count = 'Gene Count in Input, PPI, and OptDis Output'
ylabel_frac = 'Fraction of Genes'
title_frac = 'Gene Fraction in Input, PPI and OptDis Output'

# Set directories
mydirs = set_directories.mydirs(inputdir, outputdir)

def find_genes_in_input(data_fullpath, gene_symbol_colname, probe_or_gene_colname):
    '''
    From a data containing all the genes in the cohort and also which genes
    have alternative splicing events detected, count how many are alternatively spliced
    and how many total genes there are.
    Ignore NAs.
    
    Input:
        data_fullpath: containing genes fed into optdis and also whether alternative
            splice event occured.
    Output: count of how many alternative splice events occured. 
    '''
    non_as_genes = []
    as_genes = []
    NAcount = 0
    with open(data_fullpath, 'rb') as data_file:
        data_reader = csv.reader(data_file, delimiter='\t')
        header = data_reader.next()
        for row in data_reader:
            genename = row[header.index(gene_symbol_colname)]
            if genename != 'NA':
                probe_or_gene = row[header.index(probe_or_gene_colname)]
                if probe_or_gene == 'gene':
                    non_as_genes.append(genename)
                elif probe_or_gene == 'probe':
                    as_genes.append(genename)
            else:
                # print('Gene name NA found, skipping..')
                NAcount += 1
    print('No more iterations, %s NA rows skipped. ' %NAcount)
    return non_as_genes, as_genes

def find_genes_matched_to_ppi(gene_list_fullpath, non_as_genes_input, as_genes_input):
    '''
    From PPI network file (two columns of gene names), find which non_as_genes and as_genes
    match to the PPI network. 
    '''
    # Initialize output lists
    gene_list = []
    non_as_genes_ppi = []
    as_genes_ppi = []
    with open(gene_list_fullpath, 'rb') as data_file:
        data_reader = csv.reader(data_file, delimiter='\t')
        '''
        Read all genes to memory, then reduce list to unique genes only.
        See how many input genes match to ppi. 
        The list contains two columns corresponding to two gene interactions.
        '''
        for row in data_reader:
            gene_list.append(row[0])    # Gene name in column 1
            gene_list.append(row[1])    # Gene name in column 2
        gene_list = list(set(gene_list))
        for non_as_g in non_as_genes_input:
            if non_as_g in gene_list:
                non_as_genes_ppi.append(non_as_g)
        for as_g in as_genes_input:
            if as_g in gene_list:
                as_genes_ppi.append(as_g)
        
        return non_as_genes_ppi, as_genes_ppi
    
def find_genes_in_output(optdis_output_fullpath, non_as_genes_ppi, as_genes_ppi):
    '''
    From optdis results, find which genes had AS event and which ones did not. 
    Inputs:
        optdis_output_fullpath: contains an appended list of OptDis results that 
            provides whether that gene was AS or not in the 5th column.
        non_as_genes_ppi: list of genes that didnt show AS and matched to ppi.
        as_genes_ppi: genes that showed AS that matched to PPI
    Outputs:
        non_as_genes_output: list of genes were not AS in output.
        as_genes_output: list of genes showing AS in optdis output. 
    '''
    non_as_genes_output = []
    as_genes_output = []
    
    with open(optdis_output_fullpath, 'rb') as data_file:
        data_reader = csv.reader(data_file, delimiter='\t')
        header = data_reader.next()
        for row in data_reader:
            genename = row[header.index(gene_symbol_colname)]
            probe_or_gene = row[header.index(probe_or_gene_colname_optdisresults)]
            if probe_or_gene == 'gene':
                non_as_genes_output.append(genename)
            elif probe_or_gene == 'probe':
                as_genes_output.append(genename)
            else:
                sys.exit('%s must be gene or probe' %probe_or_gene)
    
    return non_as_genes_output, as_genes_output

def prepare_stacked_barplot():
    pass

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Optdis Output File and Probe/Gene Expression Data File must be specified in command line.' \
              '\nFilepath is relative to the output directory ')
        sys.exit()
    optdis_partialpath = sys.argv[1]
    probegene_partialpath = sys.argv[2]
    hprd_partialpath = sys.argv[3]
    
    optdis_fullpath = os.path.join(mydirs.outputdir, optdis_partialpath)
    probegene_fullpath = os.path.join(mydirs.outputdir, probegene_partialpath)
    hprd_fullpath = os.path.join(mydirs.outputdir, hprd_partialpath)
    
    # Classify input genes as non_AS or AS
    non_as_genes_input, as_genes_input = find_genes_in_input(probegene_fullpath, 
                                                                gene_symbol_colname, 
                                                                probe_or_gene_colname)
    # Classify ppi genes as non-AS or AS
    non_as_genes_ppi, as_genes_ppi = find_genes_matched_to_ppi(hprd_fullpath, 
                                                               non_as_genes_input, 
                                                               as_genes_input)
    # Classify optdis output genes AS non_as or AS
    non_as_genes_output, as_genes_output = find_genes_in_output(optdis_fullpath,
                                                                non_as_genes_ppi,
                                                                as_genes_ppi)
    # Plot input genes, ppi genes, optdis genes.
    non_as_genes_count = [len(non_as_genes_input), len(non_as_genes_ppi), 
                          len(non_as_genes_output)]
    as_genes_count = [len(as_genes_input), len(as_genes_ppi), 
                      len(as_genes_output)]
    print non_as_genes_count, as_genes_count
    plots.plot_stacked_bar(non_as_genes_count,
                           as_genes_count,
                           barplot_categories, barplot_events, 
                           barplot_barwidths, 
                           color_vector, xlabel, ylabel_count, title_count)
    
    # Plot input genes, ppi genes, optdis genes by proportion.
    total_count = [(i + j) for i, j in zip(non_as_genes_count, as_genes_count)]
    non_as_genes_frac = [float(i) / j for i, j in zip(non_as_genes_count, 
                                                      total_count)]
    as_genes_frac = [float(i) / j for i, j in zip(as_genes_count, total_count)]
    
    print non_as_genes_frac, as_genes_frac
    plots.plot_stacked_bar(non_as_genes_frac, as_genes_frac, 
                           barplot_categories, barplot_events, 
                           barplot_barwidths, color_vector, xlabel, 
                           ylabel_frac, title_frac)
    
    print len(non_as_genes_input), len(non_as_genes_ppi), len(non_as_genes_output), len(as_genes_input), len(as_genes_ppi), len(as_genes_output)
    
    
    
    
    
        