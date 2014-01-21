'''
Created on 2014-01-21

@author: jyeung

After t-testing and adjusting p-values,
filter events by BH-pvalue and delta psi.
'''

import sys
import csv
import os
from optparse import OptionParser

def main():
    usage = 'usage: %prog [opts] miso_t_test_file filtered_output\n'\
        'Two args must be specified in commandline: \n'\
        '1) Text file of t_test summary (bh-pvalue adjusted) of miso events\n'\
        '2) Output text file.\n'\
        'Press -h or --help for option parameter information.'
    parser = OptionParser(usage=usage)
    parser.add_option('--bh_pval_threshold', dest='bh_pval_threshold',
                      default=0.1,
                      help='BH-pval must be equal or less than '\
                        'bh_pval_threshold. Default=0.1')
    parser.add_option('--delta_psi_threshold', dest='delta_psi_threshold',
                      default=0.3,
                      help='Delta psi (abs val) must be greater or less than '\
                        'delta_psi_threshold. Default=0.3')
    parser.add_option('--bh_pval_colname', dest='bh_pval_colname',
                      default='bh_adj_pval',
                      help='Column name containing bh-pval')
    parser.add_option('--delta_psi_colname', dest='delta_psi_colname',
                      default='abs_delta_psi',
                      help='Column name containign abs value of delta psi')
    # Parse options
    (options, args) = parser.parse_args()
    
    if len(args) != 2:
        print usage
        sys.exit()
    
    miso_t_test_filepath = args[0]
    output_filepath = args[1]
    
    # Define relevant colnames of cohorts from options
    bh_pval_thres = float(options.bh_pval_threshold)
    delta_psi_thres = float(options.delta_psi_threshold)
    bh_pval_colname = options.bh_pval_colname
    delta_psi_colname = options.delta_psi_colname
    
    # Print what options are being used...
    print 'Filter options:\nBH-Pvalue <= %s\nDelta PSI >= %s' \
        %(bh_pval_thres, delta_psi_thres)
        
    # Init write file
    if os.path.isfile(output_filepath):
        print 'File %s exists.\nAborting...' \
            %output_filepath
        sys.exit()
    writefile = open(output_filepath, 'wb')
    jwriter = csv.writer(writefile, delimiter='\t')
    # read file
    writecount = 0
    with open(miso_t_test_filepath, 'rb') as readfile:
        jreader = csv.reader(readfile, delimiter='\t')
        header = jreader.next()
        # write head to file
        jwriter.writerow(header)
        for rowcount, row in enumerate(jreader):
            bh_pval = float(row[header.index(bh_pval_colname)])
            delta_psi = float(row[header.index(delta_psi_colname)])
            if bh_pval <= bh_pval_thres and delta_psi >= delta_psi_thres:
                # passes filter, write to output
                jwriter.writerow(row)
                writecount += 1
    writefile.close()
    
    print '%s rows read and %s rows written to file: %s' \
        %(rowcount, writecount, output_filepath)
    
if __name__ == '__main__':
    main()