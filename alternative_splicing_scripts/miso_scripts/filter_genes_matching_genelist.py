'''
Created on 2013-11-26

@author: jyeung

Filter a list of genes (from MISO output, e.g.) for
only genes matching a candidate set of genes.
'''

import sys
import csv
from optparse import OptionParser

def index_ipa_list(match_list):
    '''
    IPA downstream molecules file. Get gene symbols.
    '''
    gene_list = []
    genesymbol = 'Symbol'
    with open(match_list, 'rb') as readfile:
        jreader = csv.reader(readfile, delimiter='\t')
        # Skip rights reserved first line
        jreader.next()
        headers = jreader.next()
        for row in jreader:
            gene_list.append(row[headers.index(genesymbol)])
    return gene_list

def main():
    usage = 'usage: %prog miso_bf_file ar_genelist rest_genelist output_file\n'\
        'Four args must be specified in commandline: \n'\
        '1) MISO Filtered BF file.\n'\
        '2) List of AR downstream genes.\n'\
        '3) List of REST downstream genes.\n'\
        '4) Output annotated gene file.\n'
    parser = OptionParser(usage=usage)
    (_, args) = parser.parse_args()
    if len(args) < 4:
        print(usage)
        sys.exit()
    miso_bf = args[0]
    ar_downstreams = args[1]
    rest_downstreams = args[2]
    output_bf = args[3]
    
    # MISO BF gsymbol colname
    gsymbol_str = 'gsymbol'
    
    # Get lists of AR and REST
    downstream_of_ar = index_ipa_list(ar_downstreams)
    downstream_of_rest = index_ipa_list(rest_downstreams)
    
    # Init output
    outputfile = open(output_bf, 'wb')
    jwriter = csv.writer(outputfile, delimiter='\t')
    # Write headers
    jwriter.writerow(['gsymbol', 'downstream_of_AR', 'downstream_of_REST'])
    
    # Read miso bf
    with open(miso_bf, 'rb') as readfile:
        jreader = csv.reader(readfile, delimiter='\t')
        headers = jreader.next()
        ar_count = 0
        rest_count = 0
        total_count = 0
        for row in jreader:
            total_count += 1
            is_in_ar = ''
            is_in_rest = ''
            gene_symbol = row[headers.index(gsymbol_str)]
            if gene_symbol in downstream_of_ar:
                is_in_ar = 'x'
                ar_count += 1
            if gene_symbol in downstream_of_rest:
                is_in_rest = 'x'
                rest_count += 1
            jwriter.writerow([gene_symbol, is_in_ar, is_in_rest])
    print '%s in AR, %s in REST, %s total. Outputfile: %s' %(ar_count, 
                                                             rest_count, 
                                                             total_count, 
                                                             output_bf)
    
    
if __name__ == '__main__':
    main()