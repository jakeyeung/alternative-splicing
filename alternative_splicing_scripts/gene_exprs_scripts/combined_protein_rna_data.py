'''
Created on 2013-12-18

@author: jyeung
Takes LFQ data and mRNA gene exprs data
and integrates it into one text file.
'''

import sys
from optparse import OptionParser
from integrate_protein_rna_data import \
    index_lfq_data, index_mrna_data, write_lfq_mrna_data_to_file

def main():
    usage = 'usage: %prog [opt] lfq_filename gene_exprs_filename output_filename'\
        '\nThree arguments must be specified in command line:\n'\
        '1) LFQ filename, containing LFQ intensities and two replicates.\n'\
        '2) Gene exprs filename, read count data.\n'\
        '3) Output filename\n'
    parser = OptionParser(usage=usage)
    
    # colnames for lfq data
    parser.add_option('--lfq_gene_colname', dest='lfq_gene_colname',
                      default='Gene names',
                      help='Column name of gene name')
    parser.add_option('--samp1_lfq_colname1', dest='samp1_lfq_colname1',
                      default='LFQ intensity T331_1',
                      help='Column name of LFQ intensity, sample 1 replicate 1.')
    parser.add_option('--samp1_lfq_colname2', dest='sampl1_lfq_colname2',
                      default='LFQ intensity T331_2',
                      help='Column name of LFQ intensity, sample 1 replicate 2.')
    parser.add_option('--samp2_lfq_colname1', dest='sampl2_lfq_colname1',
                      default='LFQ intensity R_1',
                      help='Column name of LFQ intensity, sample 2 replicate 1.')
    parser.add_option('--samp2_lfq_colname2', dest='sampl2_lfq_colname2',
                      default='LFQ intensity R_2',
                      help='Column name of LFQ intensity, sample 2 replicate 2.')  
    # colnames for gene exprs data
    parser.add_option('--mrna_gene_colname', dest='mrna_gene_colname',
                      default='gene_name',
                      help='Column name for mRNA exprs data.')
    parser.add_option('--samp1_exprs_colname', dest='samp1_exprs_colname',
                      default='LTL331',
                      help='Column name of gene exprs for sample 1')  
    parser.add_option('--samp2_exprs_colname', dest='samp2_exprs_colname',
                      default='LTL331_R',
                      help='Column name of gene exprs for sample 2')
    (options, args) = parser.parse_args()
    
    if len(args) < 3:
        print 'Not enough args specified.\n%s' %usage
        sys.exit()
    
    lfq_filename = args[0]
    gene_exprs_filename = args[1]
    out_filename = args[2]
    
    lfq_mrna_dic = {}
    
    # Add LFQ information to dic
    lfq_mrna_dic = index_lfq_data(lfq_filename, lfq_mrna_dic, options)
    print 'lfq data indexed from file: %s' %lfq_filename
    
    # Add gene exprs to dic
    lfq_mrna_dic = index_mrna_data(gene_exprs_filename, lfq_mrna_dic, options)
    print 'mrna data indexed from file: %s' %gene_exprs_filename
    
    # Write dic to file
    write_lfq_mrna_data_to_file(lfq_mrna_dic, out_filename, options)

if __name__ == '__main__':
    main()