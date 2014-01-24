'''
Created on 2014-01-23

@author: jyeung

Reads file, gets a random subset of lines, then writes to file.
'''

import sys
from optparse import OptionParser
import random

def main():
    usage = 'usage: %prog [opts] input output\n'\
        'Two args must be specified in commandline: \n'\
        '1) Input file, first row should be column names.\n'\
        '2) Output file, a random subset of input.\n'\
        'Press -h or --help for option parameter information.'
    parser = OptionParser(usage=usage)
    parser.add_option('-n', '--n_rows', 
                      dest='n_rows',
                      help='Number of rows to extract from input. Default 548.',
                      default=548)
    (options, args) = parser.parse_args()
    
    # Convert n_rows to integer
    try:
        n_rows = int(options.n_rows)
    except ValueError:
        print '--n_rows (-n) must be an integer.'\
            '%s found.' %options.length_of_intron

    if len(args) != 2:
        print usage
        sys.exit()
    
    input_file = args[0] 
    output_file = args[1]
    
    print 'Creating subset of input: %s.\nRandom rows to extract: %s' \
        %(input_file, n_rows)
    
    with open(input_file, 'rb') as source:
        header = source.next()
        lines = [line for line in source]
        
    random_subset = random.sample(lines, n_rows)
    
    with open(output_file, 'wb') as outfile:
        # write header
        outfile.write(header)
        outfile.write(''.join(random_subset))
        
    print 'Output file: %s' %output_file
    
if __name__ == '__main__':
    main()