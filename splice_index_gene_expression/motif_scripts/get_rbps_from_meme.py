'''
Created on 2013-10-05

@author: jyeung
Parse meme output to retrieve position-specific probability matrix 
of each motif.

Returns all the RBPs that match!

Requires selenium and beautifulsoup!
'''

from optparse import OptionParser
from selenium import webdriver
from webtool_interactions.webtool_utilities import write_list_to_file
from webtool_interactions.retrieve_rbp_from_pwm import \
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

def main(meme_file, out_file, rbp_annotations):
    '''
    Parse meme file, grab probability matrix, submit to CISBP.
    '''
    # Init my objs
    meme_parser = file_parser(meme_file)
    driver = webdriver.Firefox()
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
    parser = OptionParser()
    parser.add_option('-f', '--infile', dest='meme_file', 
                      help='Txt file from meme output from which to parse.')
    parser.add_option('-r', '--rbpfile', dest='rbp_file',
                      help='Txt file from CISBP containing RBP information.')
    parser.add_option('-o', '--outfile', dest='out_file',
                      help='Filename to be saved as output.')
    parser.add_option('-a', '--annotationfile', dest='rbp_annotations',
                      help='File containing RBP annotations from CISBP database')
    (options, args) = parser.parse_args()
    main(options.meme_file, options.out_file, options.rbp_annotations)