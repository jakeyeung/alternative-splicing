'''
Created on 2013-11-14

@author: jyeung

Take a list of RBPs found in the RBPID database, then 
find the motifs of each RBP. Grab each motif and format them
so that they will be usable as input for FIMO.
'''

import sys
from optparse import OptionParser

def main():
    usage = 'usage: %prog rbp_list.txtfile rbp_db_directory '\
        'meme_db_output.outputfile\n'\
        'Three args must be specified in commandline: \n'\
        '1) List of RBPs.\n'\
        '2) RBP directory containing RBP_Information_all_motifs.txt'\
        ' and pwms_all_motifs directory.\n'\
        '3) Output file.\n'
    parser = OptionParser(usage=usage)
    (_, args) = parser.parse_args()
    if len(args) < 3:
        print 'Incorrect number of parameters specified.'
        print usage
        sys.exit()
    input_rbps_path = args[0]
    rbp_db_path = args[1]
    output_path = args[2]
    
    # Read list of RBPs, store as a list.

if __name__ == '__main__':
    main()