'''
Created on 2013-11-08

@author: jyeung
'''

import csv

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

def write_list_of_lists_to_file(llists, outputpath, header=False):
    '''
    Inputs: list of lists and output path.
    Consider each list within list of lists as a row, then write each
    row onto the output file path.
    
    If header != False, it is expected to be a list used to write header row.
    '''
    rowcount = 0
    with open(outputpath, 'wb') as writefile:
        jwriter = csv.writer(writefile, delimiter='\t')
        # Write header if it is not false.
        if header != False:
            jwriter.writerow(header)
        for l in llists:
            jwriter.writerow(l)
            rowcount += 1
    print '%s rows written to file: %s' %(rowcount, outputpath)
    return rowcount