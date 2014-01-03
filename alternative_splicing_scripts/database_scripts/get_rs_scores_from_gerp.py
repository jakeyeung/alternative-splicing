'''
Created on 2013-12-21

@author: jyeung

Read gerp score file (e.g. chr10.maf.rates)
and return RS values given genomic coordinate (chr,stop,start)
'''

import sys
import csv
import os
from optparse import OptionParser

def get_chr_list(add_chr_prefix=True):
    '''
    Create chromosome list
    [1, 2, 3, ... 22, X, Y] if add_chr_prefix==False
    [chr1, chr2, chr3, ... chr22, chrX, chrY] if add_chr_prefix==True
    All are in strings, no integers here.
    '''
    # def my prefix
    chr_str = 'chr'
    chr_list = [str(i) for i in range(1, 23)] + ['X', 'Y']
    if add_chr_prefix == True:
        chr_list = [''.join([chr_str, i]) for i in chr_list]
    return chr_list

def get_rs_scores(rates_path, start, end):
    '''
    Read rates path, return RS scores between start and end (inclusive).
    
    Inputs:
    rates_path: path to gerp scores file, which contains neutral rates 
    in first column and RS scores in second column.
    
    start: integer of genomic coordinate start
    end: integer of genomoic coordinate end.
    
    Outputs:
    rs_scores: list containing values from start to end.
    coordinates: list containing coordinates from start to end.
    '''
    # init outputs
    rs_scores = []
    coordinates = []
    
    with open(rates_path, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        for rowcount, row in enumerate(myreader):
            if rowcount < start:
                continue
            elif rowcount >= end:
                break    # no more use iterating rows.
            else:
                # get RS value, we are in the area of interest.
                rs_scores.append(row[1])
                coordinates.append(rowcount)
    return rs_scores, coordinates
            

def main():
    usage = 'usage: %prog gerp_scores_by_chr_dir chr start stop'\
        '\nFour arguments must be specified in command line:\n'\
        '1) gerp_scores_by_chr_dir: directory containing chr*.maf.rates files\n'\
        '2) Chromosome name: chr1, chr2... chr22, chrX, chrY\n'\
        '3) Genomic start coordinate (integer)\n'\
        '4) Genomic end coordinate (integer)\n'
    parser = OptionParser(usage=usage)
    (_, args) = parser.parse_args()
    if len(args) < 4:
        print 'Not enough args specified.\n%s' %usage
        sys.exit()
    
    gerp_scores_dir = args[0]
    chromo = args[1]
    try:
        start = int(args[2])
    except ValueError:
        print 'Start must be an integer. %s inputted.' %start
    try:
        end = int(args[3])
    except ValueError:
        print 'End must be an integer. %s inputted.' %end
    
    # Open the rates file, depending on the chromosome
    # First, check that user has inputed valid chromosome.
    chr_list = get_chr_list(add_chr_prefix=True)
    if chromo not in chr_list:
        print '%s not in %s.' %(chromo, chr_list)
        print 'Chromosome must be in %s.' %chr_list
        sys.exit()
    # create rates filename: chr*.maf.rates
    rates_filename = '.'.join([chromo, 'maf', 'rates'])
    # create path
    rates_path = os.path.join(gerp_scores_dir, rates_filename)
    
    # Return RS rates between start and stop
    rs_scores = get_rs_scores(rates_path, start, end)
    
    print 'RS scores between %s and %s in %s:' %(start, end, chromo)
    print rs_scores

if __name__ == '__main__':
    main()