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
from utilities import plot_utils
from scipy.stats import stats

def get_lfq_mrna_dic_keys():
    '''
    Create function to store my dic names
    This way I dont' have to hard code constanst everywhere 
    when I want to access keys in dic
    '''
    lfq_data_key = 'lfq_data'
    mrna_data_key = 'mrna_data'
    lfq_diff_key = ' lfq_diff'
    mrna_log2_fc_key = 'mrna_log2_fc'
    return lfq_data_key, mrna_data_key, lfq_diff_key, mrna_log2_fc_key

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

def index_lfq_data(lfq_filename, mydic, options, filter_out_missing_data=False):
    '''
    Read lfq filename and store as dictionary.
    dic format:
    {genename: {"lfq_intensity": lfq_diff}}
    Gene names MAY be semi-colon separated. So try to
    split and treat each gene separately.
    
    filter_out_missing_data = False by default. This means it will
    store data even if it has missing NAs.
    If True, it will only take values when samp_diff returns non-None.
    '''
    # Define column names
    samp1_lfq_colname1, \
    samp1_lfq_colname2, \
    samp2_lfq_colname1, \
    samp2_lfq_colname2, \
    gene_colname = get_lfq_colnames(options)
    
    # Define keynames
    lfq_data_key, _, lfq_diff_key, _ = \
        get_lfq_mrna_dic_keys()
    
    index_count = 0
    with open(lfq_filename, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        myheader = myreader.next()
        for row in myreader:
            # Get gene names, may be semicolon separated, so split then iterate
            genes_semicolon_sep = row[myheader.index(gene_colname)]
            genes_list = genes_semicolon_sep.split(';')
            for gene in genes_list:
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
                    mydic[gene] = {lfq_diff_key: [samp_diff], 
                                   lfq_data_key: [samp1_lfq1, 
                                                 samp1_lfq2, 
                                                 samp2_lfq1, 
                                                 samp2_lfq2]}
                    index_count += 1
                elif samp_diff == None and filter_out_missing_data == False:
                    mydic[gene] = {lfq_diff_key: ['NA'],
                                   lfq_data_key: [samp1_lfq1, 
                                                 samp1_lfq2, 
                                                 samp2_lfq1, 
                                                 samp2_lfq2]}
                else:
                    pass
    print '%s non-NA genes (from LFQ data) stored in dic' %index_count
    return mydic

def index_mrna_data(mrna_filename, mydic, options, filter_na=False):
    '''
    Read mRNA exprs data, and add gene exprs levels to dic.
    {genename: {"gene_exprs": gene_exprs_foldchange}}
    '''
    # def colnames
    samp1_exprs_colname, \
    samp2_exprs_colname, \
    gene_colname = get_mrna_colnames(options)
    
    # def keynames
    _, mrna_data_key, _, mrna_log2_fc_key = \
        get_lfq_mrna_dic_keys()
    
    rowcount = 0
    with open(mrna_filename, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        myheader = myreader.next()
        for row in myreader:
            gene = row[myheader.index(gene_colname)]
            samp1_exprs = float(row[myheader.index(samp1_exprs_colname)])
            samp2_exprs = float(row[myheader.index(samp2_exprs_colname)])
            
            # If either gene exprs are zero, go to next.
            if samp1_exprs == 0 or samp2_exprs == 0:
                if filter_na==False:
                    foldchange_log2 = 'NA'
                elif filter_na==True:
                    continue
            else:
                foldchange_log2 = \
                    math.log(float(samp2_exprs) / float(samp1_exprs), 2)
            try:
                mydic[gene].update({mrna_log2_fc_key: [foldchange_log2], 
                                    mrna_data_key: [samp1_exprs, samp2_exprs]})
            except KeyError:
                mydic[gene] = {mrna_log2_fc_key: [foldchange_log2], 
                                    mrna_data_key: [samp1_exprs, samp2_exprs]}
            rowcount += 1
    print '%s mRNA rows stored in dic.' %rowcount
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
    
    # Define keynames
    _, _, lfq_diff_key, mrna_log2_fc_key = \
        get_lfq_mrna_dic_keys()
    
    gene_colname = 'gene'
    
    raw_data = [gene_colname, samp1_lfq_colname1, samp1_lfq_colname2, 
                      samp2_lfq_colname1, samp2_lfq_colname2, 
                      samp1_exprs_colname, samp2_exprs_colname]
    processed_data = [lfq_diff_key, mrna_log2_fc_key]
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
    lfq_data_key, mrna_data_key, lfq_diff_key, mrna_log2_fc_key = \
        get_lfq_mrna_dic_keys()
        
    with open(out_fname, 'wb') as outfile:
        mywriter = csv.writer(outfile, delimiter='\t')
        # Write headers
        write_headers(mywriter, options)
        for rowcount, gene in enumerate(lfq_mrna_dic):
            try:
                lfq_list = lfq_mrna_dic[gene][lfq_data_key]
            except KeyError:
                lfq_list = 4 * ['NA']
            try:
                mrna_list = lfq_mrna_dic[gene][mrna_data_key]
            except KeyError:
                mrna_list = 2 * ['NA']
            try:
                lfq_diff = lfq_mrna_dic[gene][lfq_diff_key]
            except KeyError:
                lfq_diff = ['NA']
            try:
                mrna_diff = lfq_mrna_dic[gene][mrna_log2_fc_key]
            except KeyError:
                mrna_diff = ['NA']
            mywriter.writerow([gene] + lfq_list + mrna_list + lfq_diff + mrna_diff)
    print '%s rows written to: %s' %(rowcount, out_fname)
    return None
    
def scatter_plot_lfq_mrna(lfq_mrna_dic, spliced_genes, spliced_only=False):
    '''
    Plots log2 fold change with lfq difference. 
    '''
    # Define keynames
    _, _, lfq_diff_key, mrna_log2_fc_key = \
        get_lfq_mrna_dic_keys()
    
    # Get vectors lfq and mrna as a list.
    bubble_annotations = []
    lfq_diff_vector = []
    mrna_diff_vector = []
    color_vector = []    # blue if not spliced, red if spliced.
    for gene in lfq_mrna_dic:
        if spliced_only:
            if gene not in spliced_genes:
                continue
        if lfq_diff_key in lfq_mrna_dic[gene] and \
            mrna_log2_fc_key in lfq_mrna_dic[gene]:
            try:
                # Since these are lists of length 1, take the first element [0]
                if len(lfq_mrna_dic[gene][lfq_diff_key]) == 1:
                    lfq_diff = float(lfq_mrna_dic[gene][lfq_diff_key][0])
                else:
                    print 'Expected %s to be of length 1.'\
                        %lfq_mrna_dic[gene][lfq_diff_key]
            except ValueError:
                continue
            try:
                # Since these are lists of length 1, take the first element [0]
                if len(lfq_mrna_dic[gene][mrna_log2_fc_key]) == 1:
                        mrna_diff = \
                            float(lfq_mrna_dic[gene][mrna_log2_fc_key][0])
                else:
                    print 'Expected %s to be of length 1.'\
                        %lfq_mrna_dic[gene][mrna_log2_fc_key]
            except ValueError:
                continue
            
            if gene in spliced_genes:
                color_vector.append('r')
            else:
                # continue
                color_vector.append('b')
                
            lfq_diff_vector.append(lfq_diff)
            mrna_diff_vector.append(mrna_diff)
            bubble_annotations.append(gene)

    # plot_utils.plot_bubble_chart(mrna_diff_vector, lfq_diff_vector, bubble_annotations)
    
    # Create bubble size
    
    if spliced_only:
        plot_utils.plot_bubble_plot(mrna_diff_vector, lfq_diff_vector, color_vector, 
                         bubble_annotations, xlabel='mRNA Log 2 Fold Change', 
                         ylabel='Protein Log 2 Fold Change', 
                         title='331 vs 331R: RNA-Seq and LFQ Data', 
                         legend=['AS Genes'],
                         saveplot=False, output_fullpath=None)
    else:
        plot_utils.plot_bubble_plot(mrna_diff_vector, lfq_diff_vector, color_vector, 
                         bubble_annotations, xlabel='mRNA Log 2 Fold Change', 
                         ylabel='Protein Log 2 Fold Change', 
                         title='331 vs 331R: RNA-Seq and LFQ Data', 
                         legend=['Genes', 'AS Genes'],
                         saveplot=False, output_fullpath=None)
    
def get_spliced_genes(miso_filename):
    '''
    Get list of spliced genes from miso
    '''
    genename_list = []
    with open(miso_filename, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        myheader = myreader.next()
        for row in myreader:
            genename = row[myheader.index('gsymbol')]
            genename_list.append(genename)
    return genename_list

def split_by_splice_status(lfq_mrna_dic, spliced_genes, spliced=None):
    '''
    Given lfq mrna dic, get lfq and mrna data for three cases:
    1) Only spliced genes (spliced == True)
    2) Only non-spliced genes (spliced == False)
    3) All genes (spliced == None)
    '''
    _, _, lfq_diff_key, mrna_log2_fc_key = \
        get_lfq_mrna_dic_keys()
    
    
    lfq_diff_list = []
    mrna_diff_list = []
    
    if spliced == None:
        for gene in lfq_mrna_dic:
            if lfq_diff_key in lfq_mrna_dic[gene] and \
                mrna_log2_fc_key in lfq_mrna_dic[gene]: 
                lfq_diff_list += lfq_mrna_dic[gene][lfq_diff_key]
                mrna_diff_list += lfq_mrna_dic[gene][mrna_log2_fc_key]

    elif spliced == True:
        for gene in lfq_mrna_dic:
            if gene in spliced_genes:
                if lfq_diff_key in lfq_mrna_dic[gene] and \
                    mrna_log2_fc_key in lfq_mrna_dic[gene]:
                        lfq_diff_list += lfq_mrna_dic[gene][lfq_diff_key]
                        mrna_diff_list += lfq_mrna_dic[gene][mrna_log2_fc_key]

    elif spliced == False:
        for gene in lfq_mrna_dic:
            if gene not in spliced_genes:
                if lfq_diff_key in lfq_mrna_dic[gene] and \
                    mrna_log2_fc_key in lfq_mrna_dic[gene]: 
                        lfq_diff_list += lfq_mrna_dic[gene][lfq_diff_key]
                        mrna_diff_list += lfq_mrna_dic[gene][mrna_log2_fc_key]
            else:
                pass
    else:
        print 'Spliced argument must be None, True or False. %s found.' \
            %spliced
    return lfq_diff_list, mrna_diff_list
            
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
    parser.add_option('--spliced_only', dest='spliced_only',
                      default='False',
                      help='True or False. True shows only spliced genes. '\
                        'False shows all. Default is False.')
    (options, args) = parser.parse_args()
    
    if len(args) < 3:
        print 'Not enough args specified.\n%s' %usage
        sys.exit()
    
    lfq_filename = args[0]
    gene_exprs_filename = args[1]
    miso_filename = args[2]
    
    # parse options
    spliced_only = options.spliced_only
    if spliced_only in ['True', 'true', 'T', 'TRUE']:
        spliced_only = True
    elif spliced_only in ['False', 'false', 'F', 'FALSE']:
        spliced_only = False
    else:
        print 'Spliced only option must be True or False. %s found.' \
            %spliced_only
        sys.exit()
    print 'splicing_only: %s' %spliced_only
    
    lfq_mrna_dic = {}
    
    # Add LFQ information to dic
    lfq_mrna_dic = index_lfq_data(lfq_filename, lfq_mrna_dic, options, 
                                  filter_out_missing_data=True)
    print 'lfq data indexed from file: %s' %lfq_filename
    
    # Add gene exprs to dic
    lfq_mrna_dic = index_mrna_data(gene_exprs_filename, lfq_mrna_dic, options,
                                   filter_na=True)
    print 'mrna data indexed from file: %s' %gene_exprs_filename
    
    # Write dic to file
    # write_lfq_mrna_data_to_file(lfq_mrna_dic, out_filename, options)
    
    # Get differentially spliced genes (non-redundant only)
    spliced_genes = list(set(get_spliced_genes(miso_filename)))
    print '%s spliced genes extracted from %s' %(len(spliced_genes), 
                                                 miso_filename)
    
    # Calculate Pearson and Spearman correlation for non-AS genes and AS genes
    
    # Create x and y vectors for spliced, nonspliced and both
    spliced_mrna_log2_fc, spliced_lfq_diff = \
        split_by_splice_status(lfq_mrna_dic, spliced_genes, spliced=True)
    non_spliced_mrna_log2_fc, non_spliced_lfq_diff = \
        split_by_splice_status(lfq_mrna_dic, spliced_genes, spliced=False)
    mrna_log2_fc, lfq_diff = \
        split_by_splice_status(lfq_mrna_dic, spliced_genes, spliced=None)
    # Calculate r and pvals for Pearson

    for mrna_diff_vector, \
        lfq_diff_vector, \
        splice_status in \
            zip([spliced_mrna_log2_fc, non_spliced_mrna_log2_fc, mrna_log2_fc], 
                [spliced_lfq_diff, non_spliced_lfq_diff, lfq_diff], 
                ['DS Genes', 'Non-DS Genes', 'All Genes']):
        pearsonr, pearsonpval = \
            stats.pearsonr(mrna_diff_vector, lfq_diff_vector)
        print 'Gene set:%s\nPearson coefficient: %s\nPval:%s' \
            %(splice_status, pearsonr, pearsonpval)
        spearmanr, spearmanpval = \
            stats.spearmanr(mrna_diff_vector, lfq_diff_vector)
        print 'Gene set:%s\nSpearman coefficient: %s\nPval:%s' \
            %(splice_status, spearmanr, spearmanpval)
        slope, intercept, r_value, p_value, std_err = stats.linregress(mrna_diff_vector,lfq_diff_vector)
        print 'slope: %s\nintercept: %s\nr_value: %s\nstd_error: %s' %(slope, intercept, r_value, std_err)
        
    # Scatterplot data
    scatter_plot_lfq_mrna(lfq_mrna_dic, spliced_genes, spliced_only=spliced_only)
    
    
if __name__ == '__main__':
    main()