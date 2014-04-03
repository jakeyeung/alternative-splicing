'''
Created on 2014-03-31

@author: jyeung

Reads as input two (or more) anchor.summary files 
(created from 2 instances of run_anchor_batch.py)
and plots a density chart comparing amino acid sequence
lengths between the 2 summary files. 
'''

import sys
import os
import csv
from optparse import OptionParser
from utilities import plot_functions

def get_aa_length_from_anchor_file(jfile):
    '''
    Reads file, expects a header and looks for
    column "amino_acid_sequence" and extracts
    length of each amino acid sequence.
    
    Returns list of lengths.
    '''
    colname = 'amino_acid_sequence'
    aa_length_list = []
    with open(jfile, 'rb') as readfile:
        reader = csv.reader(readfile, delimiter='\t')
        header = reader.next()
        for row in reader:
            aa_length_list.append(len(row[header.index(colname)]))
    return aa_length_list

def main():
    usage = 'usage: %prog [opt] directory1 directory2'\
        '\nTwo arguments must be specified in command line:\n'\
        '1) Directory of interest containing .summary files\n'\
        '2) Directory to which to compare containing .summary files (control)\n'
    parser = OptionParser(usage=usage)
    (_, args) = parser.parse_args()
    
    if len(args) != 2:
        print 'Requires 2 arguments to be specified in command line'
        print usage
        sys.exit()
    # parse args
    dir0 = args[0]
    dir1 = args[1]
    
    # get all files containing *sion.summary in each directory. Those are 
    # anchor input files.
    ext = 'sion.summary'
    anchor_files0 = \
        [os.path.join(dir0, f) for f in os.listdir(dir0) if f.endswith(ext)]
    anchor_files1 = \
        [os.path.join(dir1, f) for f in os.listdir(dir1) if f.endswith(ext)]
    
    # Read each file, retrieving the lengths of each amino acid sequence
    aa_lengths0 = []
    aa_lengths1 = []
    for file0 in anchor_files0:
        aa_lengths0 += get_aa_length_from_anchor_file(file0)
    for file1 in anchor_files1:
        aa_lengths1 += get_aa_length_from_anchor_file(file1)
    
    for id, aa in zip(['xeno', 'control'], [aa_lengths0, aa_lengths1]):
        print 'Mean for %s' %id
        print sum(aa)/float(len(aa))
        print 'Median for %s' %id
        print sorted(aa)[len(aa)//2]
        print '[min,max] for %s' %id
        print '[%s,%s]' %(min(aa), max(aa))
    plot_functions.plot_density([aa_lengths0, aa_lengths1])
    
    
if __name__ == '__main__':
    main()