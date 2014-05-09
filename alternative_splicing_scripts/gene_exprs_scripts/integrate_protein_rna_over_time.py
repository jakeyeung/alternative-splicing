'''
Created on 2014-05-08

@author: jyeung

Extension of integrate_protein_rna_data.

Visualizes proteomic and mRNA data for three time points (start, middle, end)

Dotted lines connect start, middle, end for each gene, to show progression
of expressions overtime.

uses functions from integrate_protein_rna_data.py
'''

import sys
from optparse import OptionParser
import csv
from gene_exprs_scripts.integrate_protein_rna_data import is_nan

def try_float_check_nan(mystr):
    '''
    Given a string, try to convert to float unless
    it is equal to 'nan' (case insensitive)
    '''
    if mystr.lower() == 'nan':
        return mystr    # no conversion
    else:
        try:
            mystr = float(mystr)
        except ValueError:
            print 'Warning: could not float() on %s' %mystr
        return mystr
    
def get_tup_check_nan(str1, str2, nan_to_float=False):
    '''
    Given two strings, return a tuple as floats. Unless it is
    NaN, then do not convert to float. Just keep it as string.
    '''
    if nan_to_float:
        return(float(str1), float(str2))
    else:
        str1 = try_float_check_nan(str1)
        str2 = try_float_check_nan(str2)
    return (str1, str2)

def index_lfq_data_over_time(lfq_filepath,
                             colnames,
                             gene_colname,
                             nan_to_float=False):
    '''
    # Purpose:
    Returns a dictionary of form:
     
    lfq_dic = {gene: {lfq_subkey: [samp1, samp2, samp3]}...}
     
    LFQ data are in duplicates. 3 samples, 2 duplicates.
     
    Utilizes functions from integrate_protein_rna_data.py
    
    Expects colnames to be in tuples because lfq has duplicates
    example:
    colnames = [(col1a, col1b), (col2a, col2b)...]
    *
    convert_nan_to_float = False by default
    There are NaNs in LFQ data. Normally it will not convert
    these 'nan' to float unless specified.
    '''
    lfq_dic = {}
    lfq_subkey = 'lfq_intensity'
     
    with open(lfq_filepath, 'rb') as lfq_file:
        lfq_reader = csv.reader(lfq_file, delimiter='\t')
        # get headers from first row
        header = lfq_reader.next()
        for row in lfq_reader:
            # Get gene names, may be semicolon separated, so split then iterate
            genes_semicolon_sep = row[header.index(gene_colname)]
            genes_list = genes_semicolon_sep.split(';')
            for gene in genes_list:
                if gene == '':
                    # skip genes that are empty strings.
                    continue
                # init gene in subdic with empty list
                lfq_dic[gene] = {lfq_subkey: []}
                
                # get LFQ intensities for each sample. 
                # note: each sample has 2 duplicates (denoted a, b).
                # put each sample as a tuple (a, b)
                for col_pair in colnames:
                    samp_tup = \
                        get_tup_check_nan(row[header.index(col_pair[0])],
                                          row[header.index(col_pair[1])],
                                          nan_to_float=nan_to_float)
                    lfq_dic[gene][lfq_subkey].append(samp_tup)
    return lfq_dic

def index_mrna_data_over_time(mrna_filename,
                              mrna_exprs_colnames,
                              gene_colname):
    '''
    Reads mrna_filename, takes mrna expressions from
    specified columns. Outputs dictionary of form:
    
    {gene: {subkey: [samp1exprs, samp2exprs, samp3exprs...]}}
    '''
    outdic = {}
    subkey = 'mrna_intensity'
     
    with open(mrna_filename, 'rb') as mrna_file:
        mrna_reader = csv.reader(mrna_file, delimiter='\t')
        header = mrna_reader.next()
        for row in mrna_reader:
            gene = row[header.index(gene_colname)]
            # init subdic with empty list inside subkey
            outdic[gene] = {subkey: []}
            # get exprs from each sample
            for colname in mrna_exprs_colnames:
                gene_exprs = float(row[header.index(colname)])
                # put gene exprs into outdic
                outdic[gene][subkey].append(gene_exprs)
    return outdic

def count_n_nans(tup):
    '''
    Given a tuple, count how many 'NaN' are in the tuple.
    
    Also return a list of indexes locating the non-NaNs.
    
    Given tup:
    (1, 'NaN')
    
    Return:
    n_nans = 1
    number_indexes = [0]
    '''
    n_nans = 0
    number_indexes = []
    
    for i, t in enumerate(tup):
        if isinstance(t, float) and not is_nan(t):
            number_indexes.append(i)
        elif is_nan(t):
            n_nans += 1
        elif isinstance(t, str):
            if t.lower() == 'nan':
                n_nans += 1
        else:
            print 'Warning: unrecognized input: %s' %t
            print 'Expected a float, "NaN" or nan'
    return n_nans, number_indexes

def tup_to_avg(tup_lst, assign_nans='NaN'):
    '''
    LFQ data comes in duplicates for each sample.
    
    Given
    e.g.: [(samp1_dup1, samp1_dup2),..., (sampN_dup1, sampN_dup2)]
    
    Output
    [avgsamp1, avgsamp2...]
    
    where avgsamp1 = mean(samp1_dup1, samp1_dup2)
    
    this will match what mrna data looks like
    Contains options here to handle 'NaNs'
    
    If one of the two in tuple is NaN, then take value of non-NaN as mean
    
    If both are NaN, assign it as 'NaN' or a value (like a low exprs).
    
    assign_nans is option to assign to 'NaN' or some other value.
    '''
    lst = []
    for tup in tup_lst:
        # count number of NaNs in the tup and index of the NaN
        n_nans, number_i_lst = count_n_nans(tup)
        
        # three cases, 0, 1 or 2 NaNs.
        if n_nans == 0:
            # take avg
            avg = sum(tup) / float(len(tup))
        elif n_nans == 1:
            # take non-NaN value.
            avg = tup[number_i_lst[0]]
        elif n_nans == 2:
            # assign it as assign_nans
            avg = assign_nans
        else:
            print 'Warning: unknown number of NaNs. Expected 0, 1, 2'
        lst.append(avg)
    return lst

def main():
    usage = 'usage: %prog [opt] lfq_filename gene_exprs_filename output_filename'\
        '\nThree arguments must be specified in command line:\n'\
        '1) LFQ filename, containing LFQ intensities and two replicates.\n'\
        '2) Gene exprs filename, read count data.\n'
    parser = OptionParser(usage=usage)
    
    #TODO: add options later
    # colnames for lfq data
    parser.add_option('--lfq_gene_colname', dest='lfq_gene_colname',
                      default='Gene names',
                      help='Column name of gene name')
    parser.add_option('--samp1_lfq_colname1', dest='col1a',
                      default='LFQ intensity T331_1',
                      help='Column name of LFQ intensity, sample 1 replicate 1.')
    parser.add_option('--samp1_lfq_colname2', dest='col1b',
                      default='LFQ intensity T331_2',
                      help='Column name of LFQ intensity, sample 1 replicate 2.')
    parser.add_option('--samp2_lfq_colname1', dest='col2a',
                      default='LFQ intensity CCW_1',
                      help='Column name of LFQ intensity, sample 2 replicate 1.')
    parser.add_option('--samp2_lfq_colname2', dest='col2b',
                      default='LFQ intensity CCW_2',
                      help='Column name of LFQ intensity, sample 2 replicate 2.')  
    parser.add_option('--samp3_lfq_colname1', dest='col3a',
                      default='LFQ intensity R_1',
                      help='Column name of LFQ intensity, sample 2 replicate 1.')
    parser.add_option('--samp3_lfq_colname2', dest='col3b',
                      default='LFQ intensity R_2',
                      help='Column name of LFQ intensity, sample 2 replicate 2.') 
    # colnames for gene exprs data
    parser.add_option('--mrna_gene_colname', dest='mrna_gene_colname',
                      default='GeneSymbol',
                      help='Column name for mRNA exprs data.')
    parser.add_option('--samp1_exprs_colname', dest='samp1_mrna_colname',
                      default='LTL331',
                      help='Column name of gene exprs for sample 1')  
    parser.add_option('--samp2_exprs_colname', dest='samp2_mrna_colname',
                      default='LTL331_Cx_3wk',
                      help='Column name of gene exprs for sample 2')
    parser.add_option('--samp3_exprs_colname', dest='samp3_mrna_colname',
                      default='LTL331R',
                      help='Column name of gene exprs for sample 3')
    
    (options, args) = parser.parse_args()
    
    if len(args) != 2:
        print 'Not enough args specified.\n%s' %usage
        sys.exit()
    
    lfq_filename = args[0]
    gene_exprs_filename = args[1]
    
    # index lfq data
    lfq_colnames = [(options.col1a, options.col1b), 
                    (options.col2a, options.col2b), 
                    (options.col3a, options.col3b)]
    lfq_dic = index_lfq_data_over_time(lfq_filename, 
                                       lfq_colnames,
                                       options.lfq_gene_colname,
                                       nan_to_float=False)
    print 'lfq data indexed from file: %s' %lfq_filename
    
    # lfq data is in tupled list, condense it to a single value
    for gene in lfq_dic:
        # overwrite tup_lst into just a list by converting
        # tuples to a single value (the avg)
        lfq_dic[gene]['lfq_intensity'] = \
            tup_to_avg(lfq_dic[gene]['lfq_intensity'], assign_nans='NaN')
            
    # index mrna data
    gene_mrna_colnames = [options.samp1_mrna_colname, 
                          options.samp2_mrna_colname,
                          options.samp3_mrna_colname]
    mrna_dic = index_mrna_data_over_time(gene_exprs_filename,
                                         gene_mrna_colnames,
                                         options.mrna_gene_colname)
    print 'mrna data indexed from file: %s' %gene_exprs_filename
    
    # create x-y coordinates for mrna vs lfq plot (x vs y)
    for gene in lfq_dic:
        # iterate lfq_dic, not mrna_dic because there will be less iterations
        # check mrna_dic contains gene
        if gene not in mrna_dic:
            continue
        
    
if __name__ == '__main__':
    main()