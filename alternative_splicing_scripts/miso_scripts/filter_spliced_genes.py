'''
Created on 2013-12-04

@author: jyeung

Read miso file (filtered already) and I want to
keep ONLY events where the genes are sufficiently expressed
and are not hugely different (fold change let's say 4).
'''

import sys
import csv
from optparse import OptionParser
'''
# If you will plot, you will also need to import:
    import matplotlib.pyplot as plt
    import numpy as np
    import math
    from scipy.stats import gaussian_kde
'''

def index_exprs_file(exprs_file, sample_1_colname, sample_2_colname):
    '''
    Read expression file, get dic of gene name and read counts of sample1 
    and sample2.
    Needs sample1 and sample2 filenames, by convention PCa and NEPC
    '''
    gene_colname = 'gene_name'
    
    exprs_dic = {}
    with open(exprs_file, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        header = myreader.next()
        for row in myreader:
            gene = row[header.index(gene_colname)]
            samp1_readcount = row[header.index(sample_1_colname)]
            samp2_readcount = row[header.index(sample_2_colname)]
            if samp1_readcount == '0' or samp2_readcount == '0':
                # I don't want any readcounts of 0, so skip those.
                continue
            # We don't want to overwrite existing info, so check first.
            if gene not in exprs_dic:
                exprs_dic[gene] = [float(samp1_readcount), 
                                   float(samp2_readcount)]
            else:
                print '%s already exists in dictionary.' %gene
                print 'Skipping...'
    return exprs_dic

def plot_readcounts(exprs_dic):
    '''
    Exprs dic: extract floats and plot density plot.
    '''
    # Import plot specific libs
    import matplotlib.pyplot as plt
    import numpy as np
    import math
    from scipy.stats import gaussian_kde
    
    readcounts = []
    for readcount_lists in exprs_dic.values():
        # Log2 transform
        for rc in readcount_lists:
            # readcounts.append(rc)
            rc += 1
            readcounts.append(math.log(rc, 2))
        # readcounts += readcount_lists
    density = gaussian_kde(readcounts)
    xspace = np.linspace(min(readcounts), max(readcounts), 100)
    density.covariance_factor = lambda : .1
    density._compute_covariance()
    plt.plot(xspace, density(xspace))
    plt.xlabel('log2 readcounts')
    plt.ylabel('Density')
    plt.title('Distribution of read counts in log2 scale.')
    plt.show()
    return None

def convert_str_to_boolean(mystr):
    '''
    Converts str to boolean.
    '''
    if mystr in ['True', 'T', 'true', 'TRUE']:
        return True
    elif mystr in ['False', 'F', 'false', 'FALSE']:
        return False
    
def not_differentially_expressed(samp1_readcount, samp2_readcount,
                                 min_reads, max_fc):
    '''
    Check readcounts are not below min reads and
    are not above or below max fold change
    '''
    '''
    Check for two criteria:
    1) Above min reads
    2) Not large or small fold changes
    '''
    fc = samp1_readcount / samp2_readcount
    # 1) Make sure both samples are above min reads
    if samp1_readcount < min_reads or samp2_readcount < min_reads:
        # go to next row
        return False
    
    # 2) Check fold changes are not above or below threshold
    elif fc > max_fc or fc < 1.0/max_fc:
        # go to next row
        return False
    else:
        return True
    
def filter_miso_file(miso_file, exprs_dic, output_file, 
                     max_foldchange=4, 
                     min_reads=1000):
    '''
    Reads miso file, checks if exprs is above min reads and
    fold change is not greater than 4 (or the reverse, less than 1/4)
    
    If event makes it past those filters, then write to outputfile.
    '''
    # define consts
    genename_colname = 'gsymbol'
    # init outputfile
    outfile = open(output_file, 'wb')
    mywriter = csv.writer(outfile, delimiter='\t')
    
    writecount = 0
    with open(miso_file, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        header = myreader.next()
        # Write header to writefile
        mywriter.writerow(header)
        for rowcount, row in enumerate(myreader):
            genename = row[header.index(genename_colname)]
            try:
                # get readcount of samp1 and samp2
                readcount_list = exprs_dic[genename]
            except KeyError:
                # go to next row if cannot find gene in dic
                continue
            samp1_readcount = float(readcount_list[0])
            samp2_readcount = float(readcount_list[1])
            
            if not_differentially_expressed(samp1_readcount, samp2_readcount, 
                                            min_reads, max_foldchange):
                mywriter.writerow(row)
                writecount += 1
    print '%s rows written out of %s.' %(writecount, rowcount)
    print 'Output file: %s' %output_file
    
def main():
    usage = 'usage: %prog [opt] miso_file exprs_file output_miso_file\n'\
        'Three args need to be specified in command line:\n'\
        '1) miso_file after filtering.\n'\
        '2) expression file containing readcounts.\n'\
        '3) Output miso file.'
    parser = OptionParser(usage=usage)
    parser.add_option('-f', '--fold_change', dest='fold_change',
                      default=4,
                      help='Fold change before its considered too '\
                      'differentially expressed. Default is 4.')
    parser.add_option('-m', '--min_reads', dest='min_reads',
                      default=1000,
                      help='min counts before we say its too lowly exprsed '\
                      'Default is 1000.')
    parser.add_option('-s', '--sample1name', dest='sample1name',
                      default='LTL331',
                      help='Sample name of PCa sample in exprs_file.'\
                      ' Default LTL331')
    parser.add_option('-S', '--sample2name', dest='sample2name',
                      default='LTL331_R',
                      help='Sample name of NEPC sample in exprs_file.'\
                      ' Default LTL331_R')
    parser.add_option('-p', '--plot_distribution', dest='plot_distribution',
                      default='False',
                      help='True or False: plots distribution of readcounts.')
    (options, args) = parser.parse_args()
    if len(args) < 3:
        print usage
        sys.exit()
    miso_file = args[0]
    exprs_file = args[1]
    output_file = args[2]
    max_foldchange = int(options.fold_change)
    min_reads = int(options.min_reads)
    plot_distribution = convert_str_to_boolean(options.plot_distribution)
    sample1_colname = options.sample1name
    sample2_colname = options.sample2name
    
    # Index exprs file
    exprs_dic = index_exprs_file(exprs_file, sample1_colname, sample2_colname)
    
    # Plot distribution of readcounts
    if plot_distribution:
        plot_readcounts(exprs_dic)
    
    # Read miso file, decide whether to keep row or not, write to output.
    filter_miso_file(miso_file, exprs_dic, output_file, 
                     max_foldchange=max_foldchange, 
                     min_reads=min_reads)
    

if __name__ == '__main__':
    main()