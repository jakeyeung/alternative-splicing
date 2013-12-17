'''
Created on 2013-12-16

@author: jyeung

Read LFQ data and get difference between LFQ intensity 331 and 331R.
Read 331 331R gene data, get fold change (log 2 scale?)

matplotlib?
'''

import sys
import csv
import math
from optparse import OptionParser

def get_lfq_colnames(options):
    '''
    Get column names from argoptions.
    '''
    samp1_lfq_colname1 = options.samp1_lfq_colname1
    samp1_lfq_colname2 = options.sampl1_lfq_colname2
    samp2_lfq_colname1 = options.sampl2_lfq_colname1
    samp2_lfq_colname2 = options.sampl2_lfq_colname2
    genename_colname = options.lfq_gene_colname
    
    return samp1_lfq_colname1, samp1_lfq_colname2, \
            samp2_lfq_colname1, samp2_lfq_colname2, \
            genename_colname

def get_mrna_colnames(options):
    '''
    Get column name from arg options
    '''
    samp1_exprs_colname = options.samp1_exprs_colname
    samp2_exprs_colname = options.samp2_exprs_colname
    gene_colname = options.mrna_gene_colname
    return samp1_exprs_colname, samp2_exprs_colname, gene_colname
            
def is_nan(num):
    # check if it is NaN by seeing if it's equal to itself.
    return num != num
            
def try_float_get_samp_avg(replicate1, replicate2, gene):
    '''
    Given two replicates, get average of two samples.
    
    But first check if it is NaN. If it is, then return 
    the replicate that is floatable.
    
    If both are NaN, then return None.
    If one is NaN, can't do avg, jsut return the non-NaN value.
    If neither are nan, get average.
    '''
    try:
        replicate1 = float(replicate1)
        replicate2 = float(replicate2)
    except ValueError:
        pass
    
    if is_nan(replicate1) and is_nan(replicate2):
        return None
    elif not is_nan(replicate1) and is_nan(replicate2):
        return float(replicate1)
    elif is_nan(replicate1) and not is_nan(replicate2):
        return float(replicate2)
    else:
        try:
            return (float(replicate1) + float(replicate2)) / 2.0
        except ValueError:
            print 'Could not float %s or %s in gene: %s' %(replicate1, 
                                                           replicate2, 
                                                           gene)
            sys.exit()
    
def try_get_samp_diff(samp1_avg, samp2_avg):
    '''
    samp avgs could either be float or None.
    If either one is None, return None (can't do diff)
    Otherwise, return samp2_avg - samp1_avg
    '''
    if samp1_avg != None and samp2_avg != None:
        return samp2_avg - samp1_avg
    else:
        return None

def index_lfq_data(lfq_filename, mydic, options):
    '''
    Read lfq filename and store as dictionary.
    dic format:
    {genename: {"lfq_intensity": lfq_diff}}
    '''
    # Define column names
    samp1_lfq_colname1, \
    samp1_lfq_colname2, \
    samp2_lfq_colname1, \
    samp2_lfq_colname2, \
    gene_colname = get_lfq_colnames(options)
    
    with open(lfq_filename, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        myheader = myreader.next()
        for row in myreader:
            gene = row[myheader.index(gene_colname)]
            # If no gene name, go next
            if gene == '':
                continue
            '''
            For each row, extract LFQ intensities of two samples, each with
            two replicates.
            '''
            samp1_lfq1 = row[myheader.index(samp1_lfq_colname1)]
            samp1_lfq2 = row[myheader.index(samp1_lfq_colname2)]
            samp2_lfq1 = row[myheader.index(samp2_lfq_colname1)]
            samp2_lfq2 = row[myheader.index(samp2_lfq_colname2)]
            
            # get samp1 avg (or just one of the replicates, or
            #  none, depending on availability)
            samp1_avg = try_float_get_samp_avg(samp1_lfq1, samp1_lfq2, gene)
            samp2_avg = try_float_get_samp_avg(samp2_lfq1, samp2_lfq2, gene)
            
            samp_diff = try_get_samp_diff(samp1_avg, samp2_avg)
            
            # Only store floats to dic, no NaNs
            if samp_diff != None:
                mydic[gene] = {'lfq_diff': [samp_diff], 
                               'lfq_data': [samp1_lfq1, 
                                             samp1_lfq2, 
                                             samp2_lfq1, 
                                             samp2_lfq2]}
            else:
                mydic[gene] = {'lfq_diff': ['NA'],
                               'lfq_data': [samp1_lfq1, 
                                             samp1_lfq2, 
                                             samp2_lfq1, 
                                             samp2_lfq2]}
    return mydic

def index_mrna_data(mrna_filename, mydic, options):
    '''
    Read mRNA exprs data, and add gene exprs levels to dic.
    {genename: {"gene_exprs": gene_exprs_foldchange}}
    '''
    samp1_exprs_colname, \
    samp2_exprs_colname, \
    gene_colname = get_mrna_colnames(options)
    
    with open(mrna_filename, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        myheader = myreader.next()
        for row in myreader:
            gene = row[myheader.index(gene_colname)]
            samp1_exprs = float(row[myheader.index(samp1_exprs_colname)])
            samp2_exprs = float(row[myheader.index(samp2_exprs_colname)])
            
            # If either gene exprs are zero, go to next.
            if samp1_exprs == 0 or samp2_exprs == 0:
                foldchange_log2 = 'NA'
            else:
                foldchange_log2 = \
                    math.log(float(samp2_exprs) / float(samp1_exprs), 2)
            try:
                mydic[gene].update({'mrna_log2_fc': [foldchange_log2], 
                                    'mrna_data': [samp1_exprs, samp2_exprs]})
            except KeyError:
                mydic[gene] = {'mrna_log2_fc': [foldchange_log2], 
                                    'mrna_data': [samp1_exprs, samp2_exprs]}
    return mydic

def write_headers(writeobj, options):
    '''
    Write headers to writeobj
    
    Raw data means data without being processed.
    
    Processed data means calculations like foldchanges etc
    '''
    # Define column names
    samp1_lfq_colname1, \
    samp1_lfq_colname2, \
    samp2_lfq_colname1, \
    samp2_lfq_colname2, \
    _ = get_lfq_colnames(options)
    
    samp1_exprs_colname, \
    samp2_exprs_colname, \
    _ = get_mrna_colnames(options)
    
    gene_colname = 'gene'
    
    raw_data = [gene_colname, samp1_lfq_colname1, samp1_lfq_colname2, 
                      samp2_lfq_colname1, samp2_lfq_colname2, 
                      samp1_exprs_colname, samp2_exprs_colname]
    processed_data = ['lfq_diff', 'mrna_log2_fc']
    writeobj.writerow(raw_data + processed_data)
    return None

def write_lfq_mrna_data_to_file(lfq_mrna_dic, out_fname, options):
    '''
    Write dic to outfile.
    
    Make a list in order of:
    [gene samp1-lfq1, samp1-lfq2, samp2-lfq1, 
        samp2-lfq2, samp1-exprs, samp2-exprs, 
        lfq_diff, mrna_log2_fc]
    '''
    with open(out_fname, 'wb') as outfile:
        mywriter = csv.writer(outfile, delimiter='\t')
        # Write headers
        write_headers(mywriter, options)
        for rowcount, gene in enumerate(lfq_mrna_dic):
            try:
                lfq_list = lfq_mrna_dic[gene]['lfq_data']
            except KeyError:
                lfq_list = 4 * ['NA']
            try:
                mrna_list = lfq_mrna_dic[gene]['mrna_data']
            except KeyError:
                mrna_list = 2 * ['NA']
            try:
                lfq_diff = lfq_mrna_dic[gene]['lfq_diff']
            except KeyError:
                lfq_diff = ['NA']
            try:
                mrna_diff = lfq_mrna_dic[gene]['mrna_log2_fc']
            except KeyError:
                mrna_diff = ['NA']
            mywriter.writerow([gene] + lfq_list + mrna_list + lfq_diff + mrna_diff)
    print '%s rows written to: %s' %(rowcount, out_fname)
    
def scatter_plot_lfq_mrna(lfq_mrna_dic):
    '''
    Plots log2 fold change with 
    '''
    pass
            
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
    
    # Add gene exprs to dic
    lfq_mrna_dic = index_mrna_data(gene_exprs_filename, lfq_mrna_dic, options)
    
    # Write dic to file
    write_lfq_mrna_data_to_file(lfq_mrna_dic, out_filename, options)
    
    # Scatterplot data
    scatter_plot_lfq_mrna(lfq_mrna_dic)
    
if __name__ == '__main__':
    main()