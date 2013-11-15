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
from utilities import writing_utils

def read_motifs_from_file(myfile, skipheader=True):
    '''
    Given a motif file (format is A C G U followed by
    a position-weighted matrix (PWM)),
    extract the fractions as a list, then return it.
    '''
    mymotifs = []
    with open(myfile, 'rb') as readfile:
        # Skip header, contains A C G U, we know that alrdy.
        if skipheader == True:
            readfile.readline()
        for l in readfile:
            mymotifs.append(l)
    return mymotifs

def write_motif_to_file(outfile, motif, motifname, 
                         alength=4, nsites=500, evalue=0):
    '''
    After reading motif from an input file,
    write the list of motif, with motifname as header
    Default options:
    alength: number of letters, 4 since it's nucleotides.
    nsites: number of sites from which motif was enriched
    eval: E value of motif that was enriched.
    '''
    # Write motifname
    outfile.write('MOTIF %s\n' %motifname)
    # Write motif stats
    width = len(motif)
    outfile.write('letter-probability matrix: alength= %s w= %s '\
                  'nsites= %s E= %s\n' %(alength, width, nsites, evalue))
    # Write motif
    for lett_freq in motif:
        outfile.write('%s' %lett_freq)
    # Add extra space at the end...
    outfile.write('\n\n')
    return None

def main():
    usage = 'usage: %prog inputfile1 inputfile2 ... inputfileN outputfile'
    parser = OptionParser(usage=usage)
    (_, args) = parser.parse_args()
    if len(args) < 3:
        print('Two args must be specified in commandline: \n'\
              'One or more input files\n'\
              'Output file\n')
        sys.exit()
    input_file_list = args[0:-1]
    output_file = args[-1]
    
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
            motifs = read_motifs_from_file(f)
            write_motif_to_file(outfile, motifs, motifname)
            fcount += 1
    print '%s motif files crafted into: %s' %(fcount, output_file)
    
if __name__ == '__main__':
    main()