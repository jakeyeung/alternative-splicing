'''
Created on 2014-01-06

@author: jyeung
After running t_test_gene_expression.py on Beltran and VPC cohort, we now have
gene expressions and p-values stored in two textfiles (one per cohort).

Combine both text files into one (for plotting in R).

Text file will be in "long" format, one p-value and one foldchange per row
'''

import sys
import csv
from math import log
from optparse import OptionParser

def convert_csv_to_list(str_csv, convert_to_float=True):
    '''
    Convert csv to list.
    input example:
        72.94419447,54.79001181,51.51413993,9.928280808
        
    Option to convert each element to float.
    Return:
        [72, 54, 51, 9.9]
    '''
    mylist = str_csv.split(',')
    if convert_to_float:
        mylist = [float(i) for i in mylist]
    return mylist

def get_dic_keynames():
    '''
    Get keynames for dic
    '''
    cohort_str = 'cohort'
    gene_str = 'gene'
    fold_change_str = 'abs_fold_change'
    p_value_str = 'log10_p_value'
    group1_exprs_str = 'group1_exprs'
    group2_exprs_str = 'group2_exprs'
    fc_dir_str = 'fc_direction'
    return cohort_str, gene_str, fold_change_str, p_value_str, \
            group1_exprs_str, group2_exprs_str, fc_dir_str

def calculate_fold_change(group1_exprs, group2_exprs, log2=True):
    '''
    Calculate fold change between two groups.
    
    Fold change is with respect to group2, that is:
    fold change of 2 means group2 is 2 times higher than group1 for that gene
    
    If expressions are in log2 scale, log2 should be True
    If they are in linear scale, then log2 should be False.
    
    Therefore, output if log2 == True returns log2 fold change
    Output if log2 == False returns fold change on linear scale.
    '''
    group1_mean = float(sum(group1_exprs)) / len(group1_exprs)
    group2_mean = float(sum(group2_exprs)) / len(group2_exprs)
    
    if log2 == True:
        foldchange = group2_mean - group1_mean
    elif log2 == False:
        foldchange = group2_mean / group1_mean
    return foldchange

def get_pvals_fc_from_file(t_test_textfile, group1_colname, 
                           group2_colname, gene_colname, pval_colname):
    '''
    Input:
    t_test_textfile: output textfile from t_test_gene_expression.py
    
    Output:
    dictionary of form:
    {MyGene: {
    fold_change:
    p_value:}}
    '''
    # init dic keynames
    outdic = {}
    _, _, fold_change_str, p_value_str, g1_str, g2_str, _ = get_dic_keynames()
    
    with open(t_test_textfile, 'rb') as t_test_file:
        myreader = csv.reader(t_test_file, delimiter='\t')
        header = myreader.next()
        for row in myreader:
            # Get gene name
            gene_name = row[header.index(gene_colname)]
            # Get pval, as float.
            pval = float(row[header.index(pval_colname)])
            # get g1, g2, convert the CSV to a list then to float for each i.
            g1_csv = row[header.index(group1_colname)]
            g1 = convert_csv_to_list(g1_csv, convert_to_float=True)
            g2_csv = row[header.index(group2_colname)]
            g2 = convert_csv_to_list(g2_csv, convert_to_float=True)
            # Calculate fold change from g1 and g2 exprs
            fold_change = calculate_fold_change(g1, g2, log2=True)
            # Store info to dic
            if gene_name not in outdic:
                outdic[gene_name] = {}
            else:
                print 'Unexpected duplication for key: %s' %gene_name
                sys.exit()
            for key, value in \
                zip([fold_change_str, p_value_str, g1_str, g2_str], 
                    [fold_change, pval, g1_csv, g2_csv]):
                outdic[gene_name][key] = value
    return outdic

def neg_log_scale(myfloat):
    '''
    Convert float to negative log scale.
    '''
    log_trans_float = -log(myfloat, 10)
    return log_trans_float

def write_outdic_to_file(mydic, outfile, shape='long'):
    '''
    Given dic of expected form:
    
    {cohortA: {
    geneA: {
    foldchange: 10
    p_value: 0.01
    g1_exprs: 3,5,4,3
    g2_exprs: 5,9,7,1
    }}}
    
    Write to output file in long or fat format.
    By default it is long. This is useful for plotting in R
    fat is easier for reading.
    '''
    # Define keys in dic
    vpc_str = 'vpc'
    beltran_str = 'beltran'
    foldchange_key = 'foldchange'
    
    if shape == 'long':
        with open(outfile, 'wb') as writefile:
            mywriter = csv.writer(writefile, delimiter='\t')
            # Get header names
            _, gene_str, fc_str, pval_str, _, _, fc_dir_str = \
                get_dic_keynames()
            # Write header
            myheader = [gene_str, 
                        ':'.join([pval_str, vpc_str]), 
                        ':'.join([pval_str, beltran_str]), 
                        fc_str,
                        fc_dir_str]
            mywriter.writerow(myheader)
            # write dic info
            rowcount = 0
            for vpc_gene, beltran_gene in zip(mydic[vpc_str], mydic[beltran_str]):
                avg_fc = mydic[foldchange_key][fc_str]
                # We want avg_fc, absolute value
                abs_avg_fc = abs(avg_fc)
                '''
                vpc_fc = mydic[vpc_str][vpc_gene][fc_str]
                beltran_fc = mydic[beltran_str][beltran_gene][fc_str]
                avg_fc = sum([vpc_fc, beltran_fc]) / 2 
                '''
                # We want -log 10 pvalue
                vpc_pval = mydic[vpc_str][vpc_gene][pval_str]
                beltran_pval = mydic[beltran_str][beltran_gene][pval_str]
                vpc_pval_neglog = neg_log_scale(vpc_pval)
                beltran_pval_neglog = neg_log_scale(beltran_pval)
                # We want direction, +1 means g2 > g1, -1 means g1 > g2
                if avg_fc >= 0:
                    direction = 'Overexpressed'
                elif avg_fc < 0:
                    direction = 'Underexpressed'
                # Write row, must match myheader
                if vpc_gene == beltran_gene:
                    mywriter.writerow([vpc_gene, vpc_pval_neglog,
                                       beltran_pval_neglog, abs_avg_fc, 
                                       direction])
                    rowcount += 1
                else:
                    print 'Warning: %s not matching %s' %(vpc_gene, beltran_gene)
                    pass
    print '%s rows written to: %s' %(rowcount, outfile)
    return None
    
def get_fc_from_file(fold_change_filepath,
                    group1_fc_colname,
                    group2_fc_colname,
                    gene_fc_colname):
    '''
    fold_change_filepath: filepath containing expressions used for fold change
    group1_fc_colname: column name containing expressions for fold change calc1
    group2_fc_colname coulmn name containing expressions for fold change calc2
    gene_fc_colname: colum name containing gene name.
    
    Fold change is in log2 scale.
    
    Output dic of form:
    {gene: log2fc: 5}
    '''
    # def constants and dic
    _, gene_str, fold_change_str, _, _, _, _ = get_dic_keynames()
    outdic = {}
    # open file, iterate rows
    with open(fold_change_filepath, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        myheader = myreader.next()
        for row in myreader:
            gene = row[myheader.index(gene_fc_colname)]
            exprs1 = row[myheader.index(group1_fc_colname)]
            exprs2 = row[myheader.index(group2_fc_colname)]
            # calculate fold change
            fc = calculate_fold_change(exprs1, exprs2, log2=False)
            if gene not in outdic:
                outdic[gene] = {}
                outdic[gene][fold_change_str] = fc
    return outdic 
            
def main():
    usage = 'usage: %prog vpc_gene_exprs_file '\
        'beltran_t_test_file output_file\n'\
        'Requires four input arguments:\n'\
        '1) Text file of VPC cohort gene exprs\n'\
        '2) Text file of Beltran cohort gene exprs\n'\
        '3) Text file containing exprs of 2 samples/groups for fold change\n'\
        '4) Output file.\n'
    parser = OptionParser(usage=usage)
    parser.add_option('--pval_colname', dest='pval_colname',
                      default='pval',
                      help='Column name containing pvals. Default "pval".')    
    parser.add_option('--gene_colname', dest='gene_colname',
                      default='gene',
                      help='Column name containing pvals. Default "gene".')    
    parser.add_option('--group1_colname', dest='group1_colname',
                      default='group1',
                      help='Column name containing group1 exprs values. '\
                        'Default "group1".')
    parser.add_option('--group2_colname', dest='group2_colname',
                      default='group2',
                      help='Column name containing group2 exprs values. '\
                        'Default "group2".')
    (options, args) = parser.parse_args()
    
    if len(args) < 4:
        print 'Four arguments need to be specified in command line.\n'
        print usage
        sys.exit()
    vpc_filepath = args[0]
    beltran_filepath = args[1]
    fold_change_filepath = args[2]
    output_filepath = args[3]
    
    group1_colname = options.group1_colname
    group2_colname = options.group2_colname
    pval_colname = options.pval_colname
    gene_colname = options.gene_colname
    group1_fc_colname = 'LTL331'
    group2_fc_colname = 'LTL331_R'
    gene_fc_colname = 'gene_name'
    
    # store pvals to a dic.
    pval_fc_dic = {}
    for cohort, filepath in zip(['vpc', 'beltran'], 
                                [vpc_filepath, beltran_filepath]):
        pval_fc_subdic = get_pvals_fc_from_file(filepath,
                                                group1_colname,
                                                group2_colname,
                                                gene_colname, 
                                                pval_colname)
        pval_fc_dic[cohort] = pval_fc_subdic
    
    # Add fold change to dic
    pval_fc_subdic = get_fc_from_file(fold_change_filepath,
                                      group1_fc_colname,
                                      group2_fc_colname,
                                      gene_fc_colname)
    pval_fc_dic['foldchange'] = pval_fc_subdic
    

    # Write dic to file
    write_outdic_to_file(pval_fc_dic, output_filepath, shape='long')

if __name__ == '__main__':
    main()