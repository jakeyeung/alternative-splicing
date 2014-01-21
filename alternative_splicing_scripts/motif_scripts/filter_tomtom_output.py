'''
Created on 2013-12-05

@author: jyeung
Takes TomTom output file and parses it so that
it must pass certain orientation and q-value thresholds.
'''

import sys
import csv
from optparse import OptionParser

def passes_qval(row, header, qval_cutoff):
    '''
    Checks if row passes filter of min qval and strand
    '''
    # define constants
    qval_colname = 'q-value'
    
    row_qval = float(row[header.index(qval_colname)])
    if row_qval <= qval_cutoff:
        return True
    else:
        return False

def passes_orientation(row, header, filter_orientation):
    '''
    If filter_orientation is True, return true if orientation
    is +
    
    If filter_orientation is False, return true no matter what.
    Since we don't care about orientation if False.
    '''
    # def colname
    orientation_colname = 'Orientation'
    
    if filter_orientation:
        orientation = row[header.index(orientation_colname)]
        if orientation == '+':
            return True
        elif orientation == '-':
            return False
        else:
            print 'Expected + or - in orientation. %s found.' %orientation
            sys.exit()
    else:
        return True

def filter_tomtom_output(tomtom_path, output_path, qval_cutoff, 
                         filter_orientation):
    '''
    # Read TomTom file, and filter by criteria, write rows that pass
    # criteria to output path.
    '''
    # init write file
    writefile = open(output_path, 'wb')
    jwriter = csv.writer(writefile, delimiter='\t')
    writecount = 0
    
    with open(tomtom_path, 'rb') as readfile:
        jreader = csv.reader(readfile, delimiter='\t')
        # Read header, write to file
        header = jreader.next()
        jwriter.writerow(header)
        # Iterate read rows, filter rows.
        for readcount, row in enumerate(jreader):
            if passes_qval(row, header, qval_cutoff) and \
                passes_orientation(row, header, filter_orientation):
                jwriter.writerow(row)
                writecount += 1
    print '%s rows written out of %s read. Output: %s' %(writecount, 
                                                         readcount, 
                                                         output_path)
    return None

def main():
    usage = 'usage: %prog tomtom_file.txt output_file.txt\n'\
        'Two positional arguments required to '\
              'be specified in command line:\n'\
              '1) TomTom output file, containing RBP IDs.\n'\
              '2) Output file.'
    parser = OptionParser(usage=usage)
    parser.add_option('-q', '--q_cutoff', dest='qval_cutoff', default=0.15,
                      help='Specify Q-value cutoff, default 0.15.')
    parser.add_option('-s', '--filter_strand', dest='filter_strand', 
                      default='True',
                      help='Specify whether to filter strand (True searches '
                      '"+" only, False searches "+" and "-")')
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print(usage)
        sys.exit()
    tomtom_path = args[0]
    output_path = args[1]
    qval_cutoff = float(options.qval_cutoff)
    if options.filter_strand in ['True', 'true', 'TRUE', 'T', 't']:
        filter_strand = True
    elif options.filter_strand in ['False', 'false', 'FALSE', 'F', 'f']:
        filter_strand = False
    else:
        print 'Filter strand must be True or False. %s found.' \
        %options.filter_strand
        sys.exit()
    
    filter_tomtom_output(tomtom_path, output_path, 
                         qval_cutoff=qval_cutoff, 
                         filter_orientation=filter_strand)

if __name__ == '__main__':
    main()