'''
Created on 2013-10-05

@author: jyeung
Parse meme output to retrieve position-specific probability matrix 
of each motif.

Returns all the RBPs that match!

Requires selenium and beautifulsoup!
'''

import sys
import os
from optparse import OptionParser
from selenium import webdriver
from webtool_utilities import write_list_to_file
from retrieve_rbp_from_pwm import \
    query_cisrbp_get_rbp, extract_rbps_from_annotations


class file_parser(object):
    '''
    Parse a file (e.g. meme output textfile)
    '''
    def __init__(self, filepath):
        '''
        Constructor
        '''
        self.path = filepath
        self.readcount = 0
        
    def __enter__(self):
        self.file = open(self.path, 'rU')
        
    def __exit__(self, typ, value, tb):
        self.file.close()
    
    def read_next_line(self):
        '''
        Read next line in path. Requires to be in 
        with mode to work.
        '''
        self.readcount += 1
        return self.file.readline()
    
def start_of_prob_matrix(myline):
    '''
    Check if row is the header for letter-probability matrix.
    This will signal the beginning of the text we want.
    '''
    # Define constants
    header_str = 'position-specific probability matrix'    # Begin of matri
    if header_str in myline:
        return True
    else:
        return False
    
def end_of_prob_matrix(myline):
    '''
    Check if row is the END of the letter_probability matrix,
    represented by a long row of dashes.
    '''
    # Define constants
    end_str = '-------------------'
    if end_str in myline:
        return True
    else:
        return False

def main(meme_file_list, output_file_list, rbp_annotations):
    '''
    Parse meme file, grab probability matrix, submit to CISBP.
    '''
    # Init driver
    driver = webdriver.Firefox()
    
    # Begin loop
    for meme_file, out_file in zip(meme_file_list, output_file_list):
        print 'Matching motifs to CISBP for file: %s' %meme_file
        # Init my objs
        meme_parser = file_parser(meme_file)
        matched_rbps = []
        
        # Get list of known RBPs from database, useful for matching
        # with CISBP Results.
        rbp_list = extract_rbps_from_annotations(rbp_annotations)
        
        # Parse my file
        with meme_parser:
            # Iterate
            while True:
                myline = meme_parser.read_next_line()
                if myline == '':
                    print('%s rows iterated.' %meme_parser.readcount)
                    break
                if start_of_prob_matrix(myline):
                    '''
                    Row is header for matrix we want.
                    Iterate each row to collect values until
                    end of matrix.
                    '''
                    motif_name = myline
                    # Skip two lines, first line are dashes, second is header.
                    meme_parser.read_next_line()
                    meme_parser.read_next_line()
                    # Get first row in matrix.
                    myline = meme_parser.read_next_line()
                    # Iterate until end of matrix.
                    nucleotide_prob_list = []
                    while not end_of_prob_matrix(myline):
                        nucleotide_prob_list.append(myline)
                        myline = meme_parser.read_next_line()
                    # Join the list of strings into one string, adding ACGU prefix.
                    cisbp_input = ''.join([' A  C  G  U\n'] + nucleotide_prob_list)
                    print 'Matching RBPs for motif: %s' %motif_name
                    matched_rbps += query_cisrbp_get_rbp(driver, cisbp_input, rbp_list)
        matched_rbps = list(set(matched_rbps))
        # Write matched RBPs to file.
        print('RBPs found: %s.' %len(matched_rbps))
        if len(matched_rbps) > 0:
            counts = write_list_to_file(matched_rbps, out_file)
            print('%s rows written to: %s' %(counts, out_file))
    
    # Close
    driver.close()
            
if __name__ == '__main__':
    usage = 'usage: %prog [options] infile outfile rbp_annotation_file'
    parser = OptionParser(usage=usage)
    parser.add_option('-b', '--batch_mode', dest='batch_mode',
                      default=False,
                      help='Batch mode: grabs all files in same directory as '\
                      'infile and loop through each file in directory.')
    (options, args) = parser.parse_args()
    meme_file = args[0]
    out_file = args[1]
    rbp_annotations = args[2]
    
    if options.batch_mode=='True':
        # Get dir
        meme_dir = os.path.dirname(meme_file)
        # Get list of filenames in dir
        meme_filename_list = os.listdir(meme_dir)
        # Split extensions from filenames
        meme_file_noext = [os.path.splitext(i)[0] for i in meme_filename_list]
        # Add new extensions to filenames (will be our output)
        output_filename_list = [''.join([i, '.rbps']) for i in meme_file_noext]
        # Append back directory back to filenames
        meme_file_list = \
            [os.path.join(meme_dir, f) for f in meme_filename_list]
        output_file_list = \
            [os.path.join(meme_dir, o) for o in output_filename_list]
    elif options.batch_mode=='False':
        meme_file_list = [meme_file]
        output_file_list = [out_file]
    else:
        print('--batch_mode option must be True '\
              'or False. %s' %options.batch_mode)
        sys.exit()
    main(meme_file_list, output_file_list, rbp_annotations)
    
    
    
    
    
    
    
    
    
    
    
    