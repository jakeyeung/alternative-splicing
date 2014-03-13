'''
Created on 2014-02-18

@author: jyeung

Match list of genes with another list, do fisher's test
for some stats.
'''

from scipy import stats
from optparse import OptionParser
from venn_diagram_cohorts import get_events

def main():
    usage = 'usage: %prog [opts] genelist1 genelist2\n'\
        'Two args must be specified in commandline: \n'\
        '1,2) Gene lists\n'\
        'Press -h or --help for option parameter information.'
    parser = OptionParser(usage=usage)
    parser.add_option('-1', '--colname1', dest='colname1',
                      default='gsymbol',
                      help='Colname containing events in cohort1')
    parser.add_option('-2', '--colname2', dest='colname2',
                      default='gsymbol',
                      help='Colname containing events in cohort2')
    parser.add_option('--filename3', dest='filename3',
                      default=None,
                      help='Include a third list, will find intersection '\
                        'between genelist1 and genelist3 and match to '\
                            'genelist2. Default is None')
    parser.add_option('-3', '--colname3', dest='colname3',
                      default='gsymbol',
                      help='Colname containing events in cohort3. Required if '\
                        '--filename3 is not None')
    parser.add_option('-n', '--n_bg_genes', dest='n_bg',
                      default=42485,
                      help='Number of background genes. Default 42485')
    # Parse options
    (options, args) = parser.parse_args()
    
    filename1 = args[0]
    filename2 = args[1]
    
    colname1 = options.colname1
    colname2 = options.colname2
    n_bg = int(options.n_bg)
    filename3 = options.filename3
    colname3 = options.colname3
    
    gene_list1 = get_events(filename1, colname1)
    gene_list2 = get_events(filename2, colname2)
    if filename3 is not None:
        print 'Filename3 supplied: will use intersection '\
            'of\n%s\nand\n%s\n as genelist.' %(filename1, filename3)
        gene_list3 = get_events(filename3, colname3)
    
    # make everything lower case
    gene_list1 = set([gene.lower() for gene in gene_list1])
    gene_list2 = set([gene.lower() for gene in gene_list2])
    if filename3 is not None:
        gene_list3 = set([gene.lower() for gene in gene_list3])
    
    # if gene_list3 exists, find intersection with genelist1, then rename it
    # as genelist1.
    if filename3 is not None:
        gene_list1 = gene_list1 & gene_list3
    
    n_in_reg = len(gene_list2 & gene_list1)
    n_not_in_reg = len(gene_list1) - n_in_reg
    
    n_in_reg_bg = len(gene_list2) - n_in_reg
    n_not_in_reg_bg = n_bg - len(gene_list1) - n_in_reg
    
    _, pvalue = stats.fisher_exact([[n_in_reg, n_not_in_reg], [n_in_reg_bg, n_not_in_reg_bg]])
    
    print 'Out of %s genes in genelist1, %s matched to genelist2.\nPvalue: %s' \
        %(len(gene_list1), n_in_reg, pvalue)
    print 'Matched genes:'
    for g in gene_list2 & gene_list1:
        print g.upper()
        
if __name__ == '__main__':
    main()