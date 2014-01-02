'''
Created on 2013-12-31

@author: jyeung

Given chromo, start, end, retrieve string of
rs scores
'''

import sys
import csv
import os

def get_chr_list(add_chr_prefix=True):
    '''
    Create chromosome list
    [1, 2, 3, ... 22, X, Y] if add_chr_prefix==False
    [chr1, chr2, chr3, ... chr22, chrX, chrY] if add_chr_prefix==True
    All are in strings, no integers here.
    '''
    # def my prefix
    chr_str = 'chr'
    chr_list = [str(i) for i in range(23)] + ['X', 'Y']
    if add_chr_prefix == True:
        chr_list = [''.join([chr_str, i]) for i in chr_list]
    return chr_list

def open_file_get_rs(rates_path, start, end):
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
    return rs_scores
            
def get_avg_rs_score(gerp_scores_dir, chromo, start, end):
    '''
    Given gerp scores directory, chromosome start and end,
    open relevant flat text file, find chromsoome start/end and retrieve
    RS scores.
    '''
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
    
    # Get RS rates between start and stop as a list
    rs_scores = open_file_get_rs(rates_path, start, end)
    
    # Get average RS rate. 1) convert to float, 2) calc avg.
    rs_scores = [float(i) for i in rs_scores]
    rs_scores_avg = float(sum(rs_scores)) / len(rs_scores)
    
    return rs_scores_avg

