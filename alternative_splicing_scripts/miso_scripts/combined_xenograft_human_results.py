'''
Created on 2014-03-04

@author: jyeung

Quick script: it combines fasta files together
from two different results.

It makes sure no duplicate lines get copied. 
'''

import sys
import csv
import os
from optparse import OptionParser

def main():
    usage = 'usage: %prog input_directory1 input_directory2 output_directory'\
        '\nThree inputs required:\n'\
        '1) Input directory 1\n'\
        '2) Input directory 2\n'\
        '3) Output directory\n'
    parser = OptionParser(usage=usage)
    
    (_, args) = parser.parse_args()
    if len(args) != 3:
        print usage
        sys.exit()
    input_filepath1 = args[0]
    input_filepath2 = args[1]
    output_filepath = args[2]
    
    # get list of files in each path
    file_list1 = list(set(os.listdir(input_filepath1)))
    file_list2 = list(set(os.listdir(input_filepath2)))
    
    # Make sure the two lists contain same files.
    if file_list1 == file_list2:
        print 'File lists contain same files.'
    else:
        print 'Error: File lists do not contain same files.'
        print 'List 1: %s' %file_list1
        print 'List 2: %s' %file_list2
        sys.exit()
    
    for jfile in file_list1:
        writecount = 0
        # get basename for output
        base = os.path.basename(jfile)
        # Initialize output file
        outfile = \
            open(os.path.join(output_filepath, base), 'wb')
        # read first file, write everything to output
        with open(os.path.join(input_filepath1, base), 'rb') as readfile:
            '''
            # iterate two rows at a time
            # first row contains fasta header
            # second row contains sequence
            '''
            fasta_headers = []    # for appending
            for fasta_header in readfile:
                seq = readfile.next()
                fasta_headers.append(fasta_header)
                # write fastas_header and seq to file
                for line in [fasta_header, seq]:
                    outfile.write(line)
                    writecount += 1
                    
        # read second file, only write if not already written in file1
        with open(os.path.join(input_filepath2, base), 'rb') as readfile:
            '''
            # iterate two rows at a time
            # first row contains fasta header
            # second row contains sequence
            '''
            for fasta_header in readfile:
                seq = readfile.next()
                if fasta_header in fasta_headers:
                    continue    # move to next fasta_header
                else:
                    # unique fasta, write to output
                    for line in [fasta_header, seq]:
                        outfile.write(line)
                        writecount += 1
        outfile.close()
        print '%s rows written to: %s' \
            %(writecount, os.path.join(output_filepath, base))
        

if __name__ == '__main__':
    main()