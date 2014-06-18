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
    parser.add_option('-n', '--n_bg_matched', dest='n_bg_matched',
                      default=29,
                      help='Number of background genes matched to genelist2. '\
                        'Default 29 for SE.')
    parser.add_option('-N', '--n_bg_total_genes', dest='n_bg_total',
                      default=8863,
                      help='Number of background genes. Default 8863 for SE.')
    parser.add_option('-B', '--background_file', dest='bg_file',
                      default=None,
                      help='Background file of gene names. Ignores -n and -N if used.')
    parser.add_option('-b', '--background_colname', dest='bg_colname',
                      default='gsymbol',
                      help='To extract gene names from bg file, need to know '\
                        'colname corresponding to gene name. Default "gsymbol".')
    # Parse options
    (options, args) = parser.parse_args()
    
    filename1 = args[0]
    filename2 = args[1]
    
    colname1 = options.colname1
    colname2 = options.colname2
    if options.bg_file is None:
        n_bg_matched = int(options.n_bg_matched)
        n_bg_total = int(options.n_bg_total)
    else:
        bg_file = options.bg_file
        bg_colname = options.bg_colname
    filename3 = options.filename3
    colname3 = options.colname3
    
    gene_list1 = get_events(filename1, colname1)
    gene_list2 = get_events(filename2, colname2)
    if options.bg_file is not None:
        print 'Extracting genes from background file: %s' %bg_file
        bg_gene_list = get_events(bg_file, bg_colname)
        print 'Done.'
    if filename3 is not None:
        print 'Filename3 supplied: will use intersection '\
            'of\n%s\nand\n%s\n as genelist.' %(filename1, filename3)
        gene_list3 = get_events(filename3, colname3)
    
    # make everything lower case
    gene_list1 = set([gene.lower() for gene in gene_list1])
    gene_list2 = set([gene.lower() for gene in gene_list2])
    if filename3 is not None:
        gene_list3 = set([gene.lower() for gene in gene_list3])
    if options.bg_file is not None:
        bg_gene_list = set([gene.lower() for gene in bg_gene_list])
        # only consider ones not in gene_list1 because that is
        # what would happen if you were to draw without replacement.
        bg_gene_list = bg_gene_list - gene_list1
    
    # if gene_list3 exists, find intersection with genelist1, then rename it
    # as genelist1.
    if filename3 is not None:
        '''
        Uncommment this to find out your total background genes
        and background genes that match to genelist2.
        print 'Total gene set: %s' %len(gene_list1 | gene_list3)
        gene_list1 = gene_list1 | gene_list3
        '''
        gene_list1 = gene_list1 & gene_list3
    
    matched_genes = gene_list2 & gene_list1
    n_in_reg = len(matched_genes)
    n_not_in_reg = len(gene_list1) - n_in_reg
    
    if options.bg_file is None:
        n_in_reg_bg = n_bg_matched
        n_not_in_reg_bg = n_bg_total - n_bg_matched
    else:
        print 'Manually calculating matches to bg genes.'
        # manually calculate genes in background
        # since we have the gene list.
        bg_matched_genes = gene_list2 & bg_gene_list
        n_in_reg_bg = len(bg_matched_genes)
        n_not_in_reg_bg = len(bg_gene_list) - n_in_reg_bg
    
    print 'Fishers exact inputs: \nmatches to list=%s, \nnon matches '\
        'to list=%s, \nmatches to bg=%s, \nnon matches to bg=%s' \
        %(n_in_reg, n_not_in_reg, n_in_reg_bg, n_not_in_reg_bg)
    _, pvalue = stats.fisher_exact([[n_in_reg, n_not_in_reg], [n_in_reg_bg, n_not_in_reg_bg]])
    
    print 'Out of %s genes in genelist1, %s matched to genelist2.\nPvalue: %s' \
        %(len(gene_list1), n_in_reg, pvalue)
    print 'Matched genes: %s' %len(gene_list2 & gene_list1)
    for g in gene_list2 & gene_list1:
        print g.upper()
        
if __name__ == '__main__':
    main()