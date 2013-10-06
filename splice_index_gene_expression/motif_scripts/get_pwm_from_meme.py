'''
Created on 2013-10-05

@author: jyeung
Parse meme output to retrieve position-specific probability matrix 
of each motif.
'''

from optparse import OptionParser

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

def main(meme_file, out_file):
    '''
    Parse meme file, grab probability matrix, submit to CISBP.
    '''
    # Init my obj
    meme_parser = file_parser(meme_file)
    
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
                print cisbp_input
             
            
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-f', '--infile', dest='meme_file', 
                      help='Txt file from meme output from which to parse.')
    parser.add_option('-r', '--rbpfile', dest='rbp_file',
                      help='Txt file from CISBP containing RBP information.')
    parser.add_option('-o', '--outfile', dest='out_file',
                      help='Filename to be saved as output.')
    (options, args) = parser.parse_args()
    main(options.meme_file, options.out_file)