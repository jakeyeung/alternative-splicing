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

def write_meme_version(writefile, meme_version):
    '''
    Writes meme version to file.
    Adds double carriage return afterwards...
    '''
    writefile.write('MEME Version %s\n\n' %meme_version)
    return None

def write_alphabet(writefile, alphabet):
    '''
    Write alphabet parameters to file.
    Adds double carriage return.
    '''
    writefile.write('ALPHABET= %s\n\n' %alphabet)
    return None

def write_strand(writefile, strand):
    '''
    write strand and double carraige return
    '''
    writefile.write('strands: %s\n\n' %strand)
    return None

def write_bckgrd_letter_freqs(writefile, bckgrd_freqs):
    '''
    bckgrd freqs is a list, in the order of ACGT.
    Parse this and write to file.
    '''
    writefile.write('Background letter frequencies\nA %s C %s G %s T %s\n\n' \
                    %(bckgrd_freqs[0], bckgrd_freqs[1], 
                      bckgrd_freqs[2], bckgrd_freqs[3]))
    return None

def write_meme_headers(writefile, 
                       meme_version=4, 
                       alphabet='ACGT', 
                       strands='+', 
                       bckgrd_freqs=[0.25, 0.25, 0.25, 0.25]):
    '''
    Take an opened outfile (using open() function)
    and write necessary headers.
    Default options:
    meme_version = 4
    alphabet = ACGT (nucleotide)
    strands = + (can be + - if you please)
    bckgrd_freqs = frequencies of ACGT, respectively. Default 0.25 across all.
    '''
    # Write necessary information for header...
    write_meme_version(writefile, meme_version)
    write_alphabet(writefile, alphabet)
    write_strand(writefile, strands)
    write_bckgrd_letter_freqs(writefile, bckgrd_freqs)
    return None

def read_motifs_from_file(myfile):
    '''
    Given a motif file (format is A C G U followed by
    a position-weighted matrix (PWM)),
    extract the fractions as a list, then return it.
    '''
    mymotifs = []
    with open(myfile, 'rb') as readfile:
        # Skip header, contains A C G U, we know that alrdy.
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
    with open(output_file, 'wb') as outfile:
        write_meme_headers(outfile)
        # Loop input files, reading their motifs and adding it to 
        for f in input_file_list:
            motifname = os.path.basename(f)
            motifs = read_motifs_from_file(f)
            write_motif_to_file(outfile, motifs, motifname)
            
    
if __name__ == '__main__':
    main()