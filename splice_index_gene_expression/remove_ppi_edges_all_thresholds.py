'''
Created on 2013-07-04

@author: jyeung
Uses functions from remove_ppi_edges to determine number of removed edges across
a wide range of thresholds. 
'''


import sys
import os
from utilities import set_directories, plots
from remove_ppi_edges import read_gene_expression, remove_edges_from_ppi


# Set directories
mydirs = set_directories.mydirs('input', 'output')
number_of_thresholds = 20

# Set plot constants
xlabel = r'$log_2$ Gene Expression Threshold'
ylabel = 'Edges Removed'
title = 'Edges Removed vs Threshold'

# Function constants
frac_samples_to_pass_threshold = 0.75


def list_range(start, stop, length_of_list):
    '''
    From start and stop (integer) and given length of final list, 
    create a list that will go from start to stop with length_of_list
    steps. 
    '''
    list_range = [x * 1 / float(length_of_list) for x in range(0, length_of_list)]
    return list_range


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('Gene expression file, ppi_network, CSV samplenames, '\
              'fractional threshold, and output filename must be indicated in commandline.'\
              'gene expression file and ppi network must be relative to input folder.')
        sys.exit()
    gene_expression_filename = sys.argv[1]
    gene_symbol_colname = sys.argv[2]
    ppi_filename = sys.argv[3]
    samplenames_csv = sys.argv[4]
    dummy_filename = sys.argv[5]
    
    gene_expression_fullpath = os.path.join(mydirs.inputdir, gene_expression_filename)
    ppi_fullpath = os.path.join(mydirs.inputdir, ppi_filename)
    dummy_fullpath = os.path.join(mydirs.outputdir, dummy_filename)
    samplenames_list = samplenames_csv.split(',')
    
    # Find threshold by reading gene expression, and saying anything
    # below frac_threshold is too low to be considered expressed.
    gene_exprs_dic = read_gene_expression(gene_expression_fullpath, samplenames_list, gene_symbol_colname)
    
    threshold_list = list_range(0, 1, number_of_thresholds)
    frac_edges_removed = []
    gene_exprs_threshold_list = []
    
    # remove_edges_from_ppi(ppi_fullpath, gene_exprs_dic, 0.25, dummy_fullpath, header=False, write_to_file=False)
    for frac_threshold in threshold_list:
        summary_dic = remove_edges_from_ppi(ppi_fullpath, gene_exprs_dic, 
                                            frac_threshold, dummy_fullpath,
                                            frac_samples_to_pass_threshold=frac_samples_to_pass_threshold, 
                                            normalize_exprs=False,
                                            header=False, write_to_file=False)
        gene_exprs_threshold_list.append(summary_dic['gene_exprs_threshold'])
        frac_edges_removed.append(summary_dic['edges_removed'] / float(summary_dic['mapped_genes']))
    plots.plot_line_graph(gene_exprs_threshold_list, frac_edges_removed, xlabel,
                           ylabel, '%s (%s frac of samples above threshold)' \
                           %(title, frac_samples_to_pass_threshold))
    
        