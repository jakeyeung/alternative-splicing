'''
Created on 2013-11-13

@author: jyeung
Used mainly to prepare Carmen's SRRM4 motifs.
Input a number of positional arguments, then output
a MEME formatted motif file, used as input for FIMO.

Highly recommend the read input files are read-only, incase
something bad happens and it gets overwritten.
'''

import sys
import os
from optparse import OptionParser
from utilities import writing_utils, reading_utils

def main():
    usage = 'usage: %prog inputfile1 inputfile2 ... inputfileN outputfile'
    parser = OptionParser(usage=usage)
    parser.add_option('-n', '--genename', dest='genename',
                      help='Name of gene, prefix added to motifname.',
                      default='')
    (options, args) = parser.parse_args()
    if len(args) < 3:
        print('Two args must be specified in commandline: \n'\
              'One or more input files\n'\
              'Output file\n')
        sys.exit()
    input_file_list = args[0:-1]
    output_file = args[-1]
    genename = options.genename
    
    # Check output file does NOT exist
    if os.path.isfile(output_file):
        print 'Warning: output file already exists.\n'\
        'Press any key to overwrite file: %s' \
        %output_file
        raw_input()
        print 'Overwriting file: %s' %output_file
        
    # init write file
    fcount = 0
    with open(output_file, 'wb') as outfile:
        writing_utils.write_meme_headers(outfile)
        # Loop input files, reading their motifs and adding it to 
        for f in input_file_list:
            motifname = os.path.basename(f)
            # remove .txt from motifname
            motifname = ''.join(motifname.split('.')[:-1])
            motifname = ','.join([genename, motifname, 'D'])
            motifs = reading_utils.read_motifs_from_file(f)
            writing_utils.write_motif_to_file(outfile, motifs, motifname,
                                              end_line='\n\n')
            fcount += 1
    print '%s motif files crafted into: %s' %(fcount, output_file)
    
if __name__ == '__main__':
    main()