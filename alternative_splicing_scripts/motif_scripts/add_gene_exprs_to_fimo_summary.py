'''
Created on 2013-11-21

@author: jyeung
After running parse_fimo_outputs.py, we want to add two columns: one is the
P-value of the t-test for each gene, second one is whether it is OVER or UNDER
expressed in neuroendocrine prostate cancer.

We will read the file from output of t_test_gene_expression.py (preferably one
using all the genes in the CISBP-DB).
'''

import sys
import csv
from optparse import OptionParser

def index_t_test_results(t_test_output_path):
    '''
    Store t test results into a dic of form:
    {gene: {group1:[], group2:[], pval:0.2}}
    '''
    # Expected colnames:
    group1_colname = 'group1'
    group2_colname = 'group2'
    pval_colname = 'pval'
    
    # read file...
    t_test_dic = {}
    with open(t_test_output_path, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        header = myreader.next()
        for row in myreader:
            for colname in [group1_colname, group2_colname]:
                # Try to split by comma, because group1 and group2 are 
                # CSVs
                t_test_dic[colname] = row(header.index(colname.split(',')))
            # Pvalue is not a CSV, so no need to split
            t_test_dic[pval_colname] = row(header.index(pval_colname)) 
    return t_test_dic

def get_gene_exprs_direction(t_test_dic, mygene):
    '''
    Get gene expression direction from t test dic given a gene.
    Outputs either 1 or -1, where :
    1 means gene exprs is higher in 
        group 2 than in group 1.
    -1 means gene exprs is lower in
        group 2 than in group 1.
    '''
    # Expected key values (same as index_t_test_results())
    group1_colname = 'group1'
    group2_colname = 'group2'
    
    group1_gene_exprs = t_test_dic[mygene][group1_colname]
    group2_gene_exprs = t_test_dic[mygene][group2_colname]
    
    # Get average of both
    group1_gene_exprs_avg = \
        sum(group1_gene_exprs) / float(len(group1_gene_exprs))
    group2_gene_exprs_avg = \
        sum(group2_gene_exprs) / float(len(group2_gene_exprs))
    
    if group2_gene_exprs_avg - group1_gene_exprs_avg > 0:
        return 1
    elif group2_gene_exprs_avg - group1_gene_exprs_avg < 0:
        return -1
    else:
        print 'Group2 and group1 exprs averages the same for gene %s, '\
        'press any key to assign direction to be 1.' %mygene
        raw_input()
        return 1
    
def add_info_to_fimo_summary(t_test_dic, fimo_summary_path, 
                             appended_fimo_summary_path):
    '''
    Add t-test information to the parsed fimo summary path.
    Write that info to appended_fimo_summary_path
    '''
    # Expected colnames:
    rbp_name_colname = 'rbp_name'
    ''''
    # Unused colnames
    ex_intr_colname = 'exon_intron_region'
    n_occurences_colname = 'n_motif_occurences'
    med_q_val_colname = 'median_q_value'
    incl_excl_colname = 'inclusion_or_exclusion'
    '''
    
    # Expected key values (same as index_t_test_results)
    pval_colname = 'pval'    # will also be a new column name.
    
    # Define extra string called direction, will be a new colname
    dir_colname = 'direction'
    
    # Init writer
    writefile = open(appended_fimo_summary_path, 'wb')
    mywriter = csv.writer(writefile, delimiter='\t')
    
    with open(fimo_summary_path, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        header = myreader.next()
        # Add pval and direction string to header, then write to outfile.
        mywriter.writerow(header + [pval_colname, dir_colname])
        # Iterate rows, grabbing pval and direction of gene exprs difference
        for row in myreader:
            mygene = row(header.index(rbp_name_colname))
            # Get pval
            mypval = t_test_dic[mygene][pval_colname]
            # Get dictionary of differential gene expression
            gene_exprs_direction = get_gene_exprs_direction(t_test_dic, mygene)
            # Add pval and direction to row, then write to outfile.
            mywriter.writerow(row + [mypval, gene_exprs_direction])
    writefile.close()
    print 'Appended fimo summary file written to: %s'\
        %appended_fimo_summary_path

def main():
    usage = 'usage: %prog [options] t_test_output parsed_fimo_summary '\
        'appended_fimo_summary.\n'\
        'Following must be defined in command line:\n'\
        '1) T-test output containing gene names of interest '\
        '(output from t_test_gene_expression.py)\n'\
        '2) Parsed fimo summary from parse_fimo_outputs.py\n'\
        '3) Output file for appended fimo summary.\n'
    parser = OptionParser(usage=usage)
    (_, args) = parser.parse_args()
    if len(args) < 4:
        print(usage)
        sys.exit()
    t_test_output_path = args[0]
    parsed_fimo_summary_path = args[1]
    appended_fimo_summary_path = args[2]
    
    # index t-test results
    t_test_dic = index_t_test_results(t_test_output_path)
    add_info_to_fimo_summary(t_test_dic, parsed_fimo_summary_path, 
                             appended_fimo_summary_path)

if __name__ == '__main__':
    main()