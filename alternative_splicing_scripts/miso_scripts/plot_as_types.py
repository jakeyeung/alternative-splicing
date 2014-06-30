'''
Created on 2014-04-12

@author: jyeung

Given a number of filtered miso 
output files, plot the distribution of cassettes,
mutually exclusives etc...
'''

from optparse import OptionParser
import sys
import csv
from miso_scripts.utilities import plot_utils

def get_rowcount(filepath, skip_header=True, 
                 get_genenames=True, gene_colname='gsymbol'):
    '''
    Reads a file, counts how many rows are in file,
    returns counts
    '''
    genes = []
    with open(filepath, 'rb') as jfile:
        jreader = csv.reader(jfile, delimiter='\t')
        if skip_header:
            header = jreader.next()    # optional
        for rowcount, row in enumerate(jreader):
            genes.append(row[header.index('gsymbol')])
    return rowcount, genes

def main():
    usage = 'usage: %prog \n'\
        'Uses optional inputs as arguments:\n'
    parser = OptionParser(usage=usage)
    parser.add_option('-c', '--se_output', dest='se_path',
                      default=None,
                      help='MISO output of cassette exons')
    parser.add_option('-m', '--mxe_output', dest='mxe_path',
                      default=None,
                      help='MISO output of mutually exclusive exons')
    parser.add_option('-5', '--a5ss_output', dest='a5ss_path',
                      default=None,
                      help='MISO output of alt 5 prime exons')
    parser.add_option('-3', '--a3ss_output', dest='a3ss_path',
                      default=None,
                      help='MISO output of alt 3 prime exons')
    parser.add_option('-r', '--ri_output', dest='ri_path',
                      default=None,
                      help='MISO output of retained introns') 
    (options, _) = parser.parse_args()
    
    # get number of events stored into a dictionary
    counts_dic = {}
    dic_keys = ['Cassette', 'Mutually exc', "Alt 5'", "Alt 3'", 'Retained intron']
    total_counts = 0
    total_genes = []
    
    for as_type, filepath in zip(dic_keys, 
                                 [options.se_path, options.mxe_path, 
                                  options.a5ss_path, options.a3ss_path, 
                                  options.ri_path]):
        try:
            as_count, as_genes = get_rowcount(filepath, skip_header=True)
            counts_dic[as_type] = as_count
            total_counts += as_count
            total_genes += as_genes
        except TypeError:
            # maybe filepath is None because it is not specified.
            # then store as None, but alert user.
            print 'Warning: Counts of %s AS type not found.' %as_type
            print 'Unable to open file: %s' %filepath
            counts_dic[as_type] = None
    
    # get labels and yvalues for barplot
    labels = []
    as_counts = []
    colors = ['yellow', 'red', 'green', 'orange', 'blue']
    for key in dic_keys:
        labels.append(key)
        as_counts.append(counts_dic[key])
        
    plot_utils.plot_bar_plot(as_counts, labels, ylabel='Number of AS events', 
                             colors=colors,
                             title='Number of differentially '\
                                'regulated AS events (n=%s)' %total_counts, 
                             label1='AS type')
    print 'Number of AS events: %s' %total_counts
    print 'Number of genes represented: %s' %len(list(set(total_genes)))

if __name__ == '__main__':
    main()