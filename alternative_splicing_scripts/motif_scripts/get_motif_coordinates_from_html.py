'''
Created on 2013-12-22

@author: jyeung

Parse meme html output and retrieve 
start and stops of each motif hit.

Currently only works for miso events from cassette exons!
'''

import sys
import os
from optparse import OptionParser

def get_motif_start_from_motif_line(motif_line):
    '''
    Motif line example:
    chr22:19964938:19965109:-@chr22:19964229:19964246:-@chr22:19963209:19963280:-  ( 23) GGGAGGGCATGGGGG 1
    
    Split by spaces, retrieve the 3rd element from split'd line, 
    which would be "23)", then remove the character from string, ).
    
    Return 23, which should be an integer.
    '''
    motif_linesplit = motif_line.split(' ')
    motif_start = motif_linesplit[3][:-1]
    try:
        return int(motif_start)
    except ValueError:
        print 'Could not extract motif start from: %s' %motif_line
        print '%s must be an integer.' %motif_start
        sys.exit()
        
def get_motif_name_from_motif_line(motif_line):
    '''
    Motif line example:
    chr22:19964938:19965109:-@chr22:19964229:19964246:-@chr22:19963209:19963280:-  ( 23) GGGAGGGCATGGGGG 1
    
    Split by spaces, retrieve the 1st element from split line.
    Return that.
    '''
    return motif_line.split(' ')[0]

def get_seq_start_end_from_miso_event(miso_event, region_of_interest):
    '''
    miso event example:
    chr22:19964938:19965109:-@chr22:19964229:19964246:-@chr22:19963209:19963280:-
    region_of_interest: either exon_1|2|3, intron_1|2_3p|5p (cassette events only)
    
    Also need to take into account the strand + or - strand.
    '''
    pass
    return None, None

def get_region_of_interest_from_filepath(filepath):
    '''
    I expect filepath to be for example:
    ~.../fasta_unshuffled/intron_1_3p_inclusion/meme.html
    
    We want intron_1_3p to be returned.
    '''
    # Find gene region (intron or exon) based on the dir name of html file
    myregion = os.path.basename(os.path.dirname(filepath))
    
    # myregion still has _inclusion, let's remove that.
    myregion_split = myregion.split('_')
    
    # check first element in split is exon or intron
    if myregion_split[0] == 'exon' or myregion_split[0] == 'intron':
        pass
    else:
        print 'Expected exon or intron in string: %s' %myregion
        sys.exit()
    
    # rejoin, removing the last element in list
    return '_'.join(myregion_split[:-1])

def main():
    usage = 'usage: %prog meme_results_file output_file\n'\
        'Two args must be specified in commandline: \n'\
        '1) Meme HTML file containing meme results.\n'\
        '2) Output file to which results will be written.\n'
    
    parser = OptionParser(usage=usage)
    
    parser.add_option('-l', '--intron_length', dest='intron_length',
                      default=100,
                      help='Length (in basepairs) of intron sequence extracted.'\
                        'Default is 100bp from the start or end of exon.')
    
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print 'Incorrect number of parameters specified.'
        print usage
        sys.exit()
    meme_html_path = args[0]
    output_path = args[1]
    
    # def search string to know when you are reading info containing
    # where the motif was matched to the event.
    searchstring = '<input type="hidden" id="blocks'
    
    # initialize motif number, starts at 1 and increments by 1.
    motif_number = 1
    
    # append motif number to searchstring, to specify the line in which 
    # to begin extracting relevant information (miso event and start/end)
    search_motif_string = ''.join([searchstring, str(motif_number)])
    # example search motif string: "<input type="hidden" id="blocks1"
    
    # define endline
    endline = '//'    # last line signaling end of motif information
    
    # get region of interest from filename
    region_of_interest = get_region_of_interest_from_filepath(meme_html_path)
    
    # Read meme html path, retrieve genomic coordinates of the motif region
    # Create search string to match 
    with open(meme_html_path, 'rb') as readfile:
        for line in readfile:
            if not line.startswith(search_motif_string):
                continue    # not the line we're looking for, keep going.
            else:
                '''
                # This is the line we're looking for, so skip one more line
                # (next line is BL    MOTIF X width=w seqs=n) and begin
                extracting miso event and start from the line.
                '''
                # skip one more line, and begin extraction
                readfile.next()
                motif_line = readfile.next()
                while not motif_line.startswith(endline):
                    miso_event = get_motif_name_from_motif_line(motif_line)
                    motif_start = get_motif_start_from_motif_line(motif_line)
                    
                    seq_start, seq_end = \
                        get_seq_start_end_from_miso_event(miso_event, 
                                                          region_of_interest)
                    
                    print miso_event
                    print motif_start
                    
                    motif_line = readfile.next()
                    raw_input()
            
    
    
if __name__ == '__main__':
    main()