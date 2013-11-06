'''
Created on 2013-07-05

@author: jyeung

Plot histogram of gene expression to get a sense of how its distribution.
'''


import os
import csv
import sys
from utilities import set_directories, plots


# Set directories
mydirs = set_directories.mydirs('input', 'output')

# Set plot parameters
xlabel = r'$log_2\ Gene\ Expression$'
ylabel = r'$Frequency$'


def get_gene_exprs_from_file(gene_expression_fullpath, samplenames_list):
    '''
    Read gene expression file, create list of gene expressions using data values
    corresponding to samplenames_list only (so a subset or all of samples).
    
    The gene expression file requires headers with samplenames, which will be used
    to locate the column index to record relevant gene expressions.
    '''
    with open(gene_expression_fullpath, 'rb') as readfile:
        reader = csv.reader(readfile, delimiter='\t')
        # Get headers as first row/
        colnames = reader.next()
        # Now iterate rest of rows, getting samples's gene exprs.
        gene_exprs_list = []
        for row in reader:
            for s in samplenames_list:
                gene_exprs_list.append(float(row[colnames.index(s)]))
    return gene_exprs_list    


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Gene expression file, CSV samplenames, '\
              'and plot title must be indicated in commandline.'\
              'gene expression file and ppi network must be relative to input folder.'\
              '\nPlot title words should be separated by underscore.')
        sys.exit()
    gene_expression_filename = sys.argv[1]
    samplenames_csv = sys.argv[2]
    plot_title_underscore = sys.argv[3]
    
    # Replace underscores with spaces
    plot_title = plot_title_underscore.replace('_', ' ')
    
    gene_expression_fullpath = os.path.join(mydirs.inputdir, gene_expression_filename)
    # plot_fullpath = os.path.join(mydirs.outputdir, plot_filename)
    samplenames_list = samplenames_csv.split(',')
    
    gene_exprs_list = get_gene_exprs_from_file(gene_expression_fullpath, samplenames_list)
    
    plots.plot_histogram(gene_exprs_list, xlabel, ylabel, plot_title, nbins=60)
    
    