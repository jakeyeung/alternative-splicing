'''
Created on 2013-10-20

@author: jyeung

Calculate differentially expressed RBPs or any list of genes.
'''

from optparse import OptionParser
import sys
import csv
from scipy import stats

def extract_column_from_table(myfile, column_index_to_extract=1, header=False):
    '''
    Grab gene names from myfile, which is likely
    an output from motif_scripts/get_rbps_from_table.py
    
    Extracts from second column, by default.
    Assumes no header in filename, by default.
    '''
    mylist = []
    with open(myfile, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        if header==True:
            myreader.next()    # Skip header if there is header.
        for row in myreader:
            mylist.append(row[column_index_to_extract])
    return mylist

def index_gene_exprs(gene_exprs_fname, gene_list, group_1, group_2, 
                     gene_colname):
    '''
    Read gene expressions, assumes first column is gene names,
    other columns contain sample names found in group_1 and group_2.
    We want a dictionary with the form:
    {gene: {group1: [samp1, samp2...], group2: [sampA, sampB...]}}
    
    If dealing with Rubin cohort...
    Group 1 should be PCa samples, Group 2 should be NEPC...
    '''
    # Initialize column names
    # gene_colname = 'gene'
    passcount = 0
    
    # Initialize keynames as gene names.
    gene_exprs_dic = {}
    for g in gene_list:
        gene_exprs_dic[g] = {'group1':[], 'group2':[]}
    
    # Initialize read object for gene exprs..
    with open(gene_exprs_fname, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        header = myreader.next()    # has a header.
        for row in myreader:
            row_gene = row[header.index(gene_colname)]
            # Put gene expression into respective group 1 and group 2.
            # Be warned, if you specified gene list, then row_gene
            # may not be in your dictionary. So check first before
            # putting gene expression into groups.
            if row_gene in gene_exprs_dic:
                for samp in group_1:
                    gene_exprs_dic[row_gene]['group1'].\
                        append(row[header.index(samp)])
                for samp in group_2:
                    gene_exprs_dic[row_gene]['group2'].\
                        append(row[header.index(samp)])
            else:
                passcount += 1
    print '%s genes not considered in this t-test.' %passcount
    return gene_exprs_dic

def t_test_gene_exprs_dic(gene_exprs_dic):
    '''
    The dictionary contains the form:
    {gene: {group1: [samp1, samp2...], group2: [sampA, sampB...]}}
    We want to loop each gene, and t-test the two groups, the p-value 
    will be appended as a new subkey.
    
    Example Output dic:
    {gene: {group1: [samp1, samp2...], group2: [sampA, sampB...], pval: 0.05}}
    '''
    emptycount = 0
    for gene in gene_exprs_dic.keys():
        # Only t_test non-empty lists...
        if not len(gene_exprs_dic[gene]['group1']) == 0 |\
                    len(gene_exprs_dic[gene]['group2']) == 0:
            g1_exprs = [float(i) for i in gene_exprs_dic[gene]['group1']]
            g2_exprs = [float(i) for i in gene_exprs_dic[gene]['group2']]
            _, pval = stats.ttest_ind(g1_exprs, g2_exprs)
            gene_exprs_dic[gene]['pval'] = [pval]
        else:
            # Prevents errors in write_dic_to_file by creating empty pval subkey 
            gene_exprs_dic[gene]['pval'] = []    
            emptycount += 1
    print '%s genes did not have gene exprs info.' %emptycount  
    return gene_exprs_dic

def write_dic_to_file(mydic, mykeys, mysubkeys, output_file):
    '''
    Writes dic of the form: 
    {gene: {group1: [samp1, samp2...], group2: [sampA, sampB...], pval: 0.05}}
    to file.
    
    Inputs:
        mydic
        mykeys such as ['gene']
        mysubkeys such as ['group1', 'group2', 'pval']
    
    Outputs:
        File containing keys and subkeys as header,
        then gene, g1_samp_exprs, g2_samp_exprs, pval
        With lists condensed to CSV format.
    '''
    # Init write counts
    writecounts = 0
    
    # Initialize output file
    with open(output_file, 'wb') as writefile:
        mywriter = csv.writer(writefile, delimiter='\t')
        # Write Header
        mywriter.writerow(mykeys + mysubkeys)
    
        # Loop through each gene, writing each gene info as a row.
        for k in mydic.keys():
            # Each gene is one row, so initialize empty list for each gene
            row = []    # Append it with subkey info.
            row.append(k)
            for subkey in mysubkeys:
                '''
                # Condense to CSV, it's ok even for pval because
                # pval was set to be a list of length 1, which 
                # will have no effect during ',' joining.
                '''
                # Convert to string before making CSV...
                list_as_string = [str(i) for i in mydic[k][subkey]]
                row.append(','.join(list_as_string))
            mywriter.writerow(row)
            writecounts += 1
    return writecounts

def main():
    usage = 'usage: %prog [optional_gene_list_file] gene_exprs_file '\
        'group_1_sampnames_file group_2_sampenames_file output_file'
    parser = OptionParser(usage=usage)
    parser.add_option('-l', '--gene_list', dest='gene_list_file',
                      default=None,
                      help='Gene list, no headers.')
    parser.add_option('-c', '--col_index', dest='col_index',
                      default=1,
                      help='Column index where gene names are. 0 is first column.')
    parser.add_option('-g', '--gene_colname', dest='gene_colname',
                      default='gene',
                      help='Column name containing gene names in exprs file.\n'\
                        'Default "gene"')    
    (opts, args) = parser.parse_args()
    
    if len(args) < 4:
        print('Following must be specified in command line:\nGene Expression '\
              'Filename\nGroup 1 Samples Filename\nGroup 2 Samples '\
              'Filename\nOutput Filename\n-h for help.')
        sys.exit()
    # args
    gene_exprs_fname = args[0]
    group_1_samps_fname = args[1]
    group_2_samps_fname = args[2]
    output_fname = args[3]
    # options
    gene_list_fname = opts.gene_list_file
    gene_colname = opts.gene_colname
    if gene_list_fname != None:
        col_index = int(opts.col_index)
    
    # If gene_list is specified, use that to get gene_list, otherwise
    # use the entire gene_exprs_fname to get gene_list.
    if gene_list_fname != None:
        gene_list = extract_column_from_table(gene_list_fname,
                                        column_index_to_extract=col_index,
                                        header=False)
    else:
        gene_list = extract_column_from_table(gene_exprs_fname, 
                                        column_index_to_extract=0, 
                                        header=True)
    print '%s genes to be tested for differential expression' %len(gene_list)
    
    # Grab group 1 and group 2 sample names as a list.
    group_1 = extract_column_from_table(group_1_samps_fname, 
                                        column_index_to_extract=0, 
                                        header=False)
    group_2 = extract_column_from_table(group_2_samps_fname, 
                                        column_index_to_extract=0, 
                                        header=False)
    # Index gene expression values of group1 and group2, all as a dictionary.
    gene_exprs_dic = index_gene_exprs(gene_exprs_fname, 
                                      gene_list,
                                      group_1,
                                      group_2,
                                      gene_colname)
    
    # Do a t-test for each gene, append it to gene_exprs_dic
    gene_exprs_dic = t_test_gene_exprs_dic(gene_exprs_dic)
    
    # Write dic to file
    writecounts = write_dic_to_file(gene_exprs_dic, ['gene'],
                                    ['group1', 'group2', 'pval'], 
                                    output_fname)
    print '%s rows written to: %s' %(writecounts, output_fname)
    
if __name__ == '__main__':
    main()