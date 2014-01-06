'''
Created on 2014-01-06

@author: jyeung
After running t_test_gene_expression.py on Beltran and VPC cohort, we now have
gene expressions and p-values stored in two textfiles (one per cohort).

Combine both text files into one (for plotting in R).

Text file will be in "long" format, one p-value and one foldchange per row
'''

import sys
from optparse import OptionParser

def get_t_test_dic(t_test_textfile):
    '''
    Input:
    t_test_textfile: output textfile from t_test_gene_expression.py
    
    Output:
    dictionary of form:
    {gene: 
    fold_change:
    p_value:}
    '''
    pass

def main():
    usage = 'usage: %prog vpc_gene_exprs_file '\
        'beltran_t_test_file output_file\n'\
        'Requires four input arguments:\n'\
        '1) Text file of VPC cohort gene exprs\n'\
        '2) Text file of Beltran cohort gene exprs\n'\
        '3) Output file.\n'
    parser = OptionParser(usage=usage)    
    (options, args) = parser.parse_args()
    
    if len(args) < 3:
        print 'Three arguments need to be specified in command line.\n'
        print usage
        sys.exit()
    vpc_filepath = args[0]
    beltran_filepath = args[1]
    output_filepath = args[2]
    
    # get t_test dic for VPC and Beltran
    t_test_dic = get_t_test_dic()

if __name__ == '__main__':
    main()