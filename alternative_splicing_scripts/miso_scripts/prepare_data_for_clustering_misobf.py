'''
Created on 2013-11-25

@author: jyeung
Parses filtered BF comparison between two samples
and returns a matrix of PSI values for both samples, 
ready to be heatplotted.
'''

import sys
import csv
from optparse import OptionParser


def main():
    # Get parse options.
    usage = 'usage: %prog [options] miso_bf_results write_output_file.\n'\
        'Following must be specified in command line:\n'\
              'miso_bf_results: output from filtered PSI comparison'\
                ' from two samples.\n'\
              'write_output_file.\n'
    parser = OptionParser(usage=usage)
    parser.add_option('-S', '--sample1_name', dest='sample1_name',
                      help='Name of sample1, default is sample1.',
                      default='sample1')
    parser.add_option('-s', '--sample2_name', dest='sample2_name',
                      help='Name of sample2, default is sample2.',
                      default='sample2')
    (opts, args) = parser.parse_args()
    if len(args) < 2:
        print(usage)
        sys.exit()
    # Define args to variables
    bf_file = args[0]
    output_file = args[1]
    sample1_name = opts.sample1_name
    sample2_name = opts.sample2_name
    
    # Define read column names
    gsymbol_str = 'gsymbol'    # also a write column name.
    samp1_psi_med_str = 'sample1_posterior_mean'
    samp2_psi_med_str = 'sample2_posterior_mean'
    
    # Prepare write file with headers.
    writecount = 0
    with open(output_file, 'wb') as writefile:
        mywriter = csv.writer(writefile, delimiter='\t')
        # Write header
        mywriter.writerow([gsymbol_str, sample1_name, sample2_name])
        # Read file, then iterate through rows
        with open(bf_file, 'rb') as readfile:
            myreader = csv.reader(readfile, delimiter='\t')
            # get headers from reader
            readheader = myreader.next()
            for row in myreader:
                genename = row[readheader.index(gsymbol_str)]
                samp1_med_psi = row[readheader.index(samp1_psi_med_str)]
                samp2_med_psi = row[readheader.index(samp2_psi_med_str)]
                # Write info to write file
                mywriter.writerow([genename, samp1_med_psi, samp2_med_psi])
                writecount += 1
    print '%s rows written to file: %s' %(writecount, output_file)

if __name__ == '__main__':
    main()