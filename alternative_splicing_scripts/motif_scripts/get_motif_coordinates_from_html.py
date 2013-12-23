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

def get_strand_from_miso_event(miso_event):
    '''
    Given event: 
    chr22:19964938:19965109:-@chr22:19964229:19964246:-@chr22:19963209:19963280:-
    Return either '-' or '+'
    '''
    return miso_event.split(':')[-1]

def get_chr_from_miso_event(miso_event):
    '''
    Given event:
    chr22:19964938:19965109:-@chr22:19964229:19964246:-@chr22:19963209:19963280:-
    Return chr22
    '''
    return miso_event.split(':')[0]

def get_intron_starts_ends_list(miso_event):
    '''
    Given event:
    chr22:19964938:19965109:-@chr22:19964229:19964246:-@chr22:19963209:19963280:-
    Returns all possible intron starts and ends
    in this caseA:
    starts = [19964246, 19965109]
    ends = [19963209, 19964229]
    '''
    strand = get_strand_from_miso_event(miso_event)
    intron_starts_list = []
    intron_ends_list = []
    # split by @ to get exon starts and stops
    miso_at_split = miso_event.split('@')
    if strand == '+':
        '''
        # iterate from i to length - 1, because
        the last exon does not have an intron start, only an end.
        We match start with exon_end in i position
        We match end with exon_start in i+1 position
        Example:
        miso_at_split:
        ['chr2:131890475:131890564:+', 'chr2:131897740:131897848:+', 'chr2:131904210:131907425:+']
        Intron starts are:
        131890564, 131897740
        Intron ends are:
        131897848, 131904210
        '''
        for i in xrange(0, len(miso_at_split) - 1):
            current_exon_coords = miso_at_split[i].split(':')
            next_exon_coords = miso_at_split[i + 1].split(':')
            intron_starts_list.append(current_exon_coords[2])
            intron_ends_list.append(next_exon_coords[1])
            
    elif strand == '-':
        '''
        # iterate from i to length - 1, because
        the last exon does not have an intron end, only an start.
        We match start with exon_end in i+1 position
        We match end with exon_start in i position
        Example:
        miso_at_split:
        ['chr22:19964938:19965109:-', 'chr22:19964229:19964246:-', 'chr22:19963209:19963280:-']
        Intron starts are:
        19964246, 19963280
        Intron ends are:
        19964938, 19964229
        '''   
        for i in xrange(0, len(miso_at_split) - 1):
            current_exon_coords = miso_at_split[i].split(':')
            next_exon_coords = miso_at_split[i + 1].split(':')
            intron_starts_list.append(next_exon_coords[2])
            intron_ends_list.append(current_exon_coords[1])
    
    else:
        print 'Expected strand to be "+" or "-". %s found.' %strand
        sys.exit()
    
    '''
    # Convert starts ends to int.
    # Starts must be incremented by ONE because technically the start
    # is an exon end, so we want the basepair after exon end.
    # Ends must be deincremented by ONE so it does not overlap
    # with exon start
    '''
    try:
        intron_starts = [int(i) + 1 for i in intron_starts_list]
    except ValueError:
        print 'Could not convert elements in %s to int.' %intron_starts_list
    try:
        intron_ends = [int(i) - 1 for i in intron_ends_list]
    except ValueError:
        print 'Could not convert elements in %s to int.' %intron_ends_list
    return intron_starts, intron_ends

def get_intron_starts_ends(miso_event, 
                           intron_number,
                           intron_3p_or_5p,
                           seq_length):
    '''
    Given event, retrieve intron starts and ends.
    Example:
    chr22:19964938:19965109:-@chr22:19964229:19964246:-@chr22:19963209:19963280:-
    if intron_number = 1, intron_3p_or_5p = "3p", seq_length = 100
    Return:
    start = 19964838 (19964938 minus 100)
    end = 19964938 
    
    Example2:
    chr2:131890475:131890564:+@chr2:131897740:131897848:+@chr2:131904210:131907425:+
    if intron_number = 2, intron_3p_or_5p = "5p", seq_length = 50
    Return:
    start = 131897848
    end = 131897898 (131897848 + 50)
    
    Intron number should be int
    '''
    strand = get_strand_from_miso_event(miso_event)
    intron_starts, intron_ends = get_intron_starts_ends_list(miso_event)
    
    # get seq length - 1 because start coordinates are inclusive.
    seq_length_adj = seq_length - 1
    
    # get intron start and ends corresponding to intron index
    intron_index = intron_number - 1

    # Get seq length at 3p or 5p site of intron start/end
    if intron_3p_or_5p == '5p':
        intron_start = intron_starts[intron_index]
        intron_end = intron_start + seq_length_adj
    elif intron_3p_or_5p == '3p':
        intron_end = intron_ends[intron_index]
        intron_start = intron_end - seq_length_adj
    else:
        print 'Expected intron_3p_or_5p to be "5p" or "3p". %s found.' \
            %intron_3p_or_5p
        sys.exit()
    return intron_start, intron_end
        
def get_exon_starts_ends(miso_event, exon_number, seq_length):
    '''
    #TODO
    '''
    pass
    
def get_seq_start_end_from_miso_event(miso_event, region_of_interest, 
                                      seq_lengths_dic):
    '''
    miso event example:
    chr22:19964938:19965109:-@chr22:19964229:19964246:-@chr22:19963209:19963280:-
    region_of_interest: either exon_1|2|3, intron_1|2_3p|5p (cassette events only)
    
    Also need to take into account the strand + or - strand.
    
    Code each region of interest separately
    possible regions of interest: exon_1|2|3, intron_1|2_3p|5p
    
    Return chr, start, end
    '''
    region_split = region_of_interest.split('_')
    exon_or_intron = region_split[0]
    
    strand = get_strand_from_miso_event(miso_event)
    chromo = get_chr_from_miso_event(miso_event)
    seq_length = seq_lengths_dic[miso_event]
    
    # check if exon or intron
    if exon_or_intron == 'exon':
        start, end = get_exon_starts_ends(miso_event, strand)
    elif exon_or_intron == 'intron':
        # check which intron number and 3p or 5p
        try:
            intron_number = int(region_split[1])
        except ValueError:
            print 'Expected int for %s' %intron_number
            sys.exit()
        intron_3p_or_5p = region_split[2]
        start, end = get_intron_starts_ends(miso_event, 
                                            intron_number, 
                                            intron_3p_or_5p,
                                            seq_length)
        print miso_event, intron_number, intron_3p_or_5p, seq_length
        print start, end
        raw_input()
        raw_input()
    else:
        print 'Expected %s to begin with "exon" or "intron". %s found.' %(region_of_interest, exon_or_intron)
        sys.exit()
    return chromo, None, None

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

def get_seq_lengths_from_meme_file(meme_file):
    '''
    Read meme file, and find the text block where it contains
    information of miso event
    '''
    # define string that signals the beginning and end relevant textblock
    begin_string = '<input type="hidden" name="combinedblock" value="'
    end_string = '">'    # signals end of textblock
    
    seq_lengths_dic = {}
    
    with open(meme_file, 'rb') as readfile:
        for line in readfile:
            if not line.startswith(begin_string):
                continue
            else:
                relevant_line = readfile.next()
                while not relevant_line.startswith(end_string):
                    '''
                    Relevant line looks like:
                    chr2:25650404:25650500:-@chr2:25642384:25642404:-@chr2:25611071:25611230:- 6.54e-04 1 100 +3 74 2.30e-05 
                    Split by ' ', retrieve "100" (4th element in split'd list)
                    '''
                    line_split = relevant_line.split(' ')
                    miso_event = line_split[0]
                    seq_length = line_split[3]
                    if miso_event not in seq_lengths_dic:
                        try:
                            seq_lengths_dic[miso_event] = int(seq_length)
                        except ValueError:
                            print 'Expected %s to be integer.' %seq_length
                            sys.exit()
                    else:
                        print 'Unexpected duplicate of %s.' %miso_event
                        print 'Press enter to overwrite.'
                        raw_input()
                    relevant_line = readfile.next()
    return seq_lengths_dic

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
    
    # create dictionary of miso_event: sequence length from meme html file
    seq_lengths_dic = get_seq_lengths_from_meme_file(meme_html_path)

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
                    
                    chromo, seq_start, seq_end = \
                        get_seq_start_end_from_miso_event(miso_event, 
                                                          region_of_interest,
                                                          seq_lengths_dic)                    
                    motif_line = readfile.next()
    
    
if __name__ == '__main__':
    main()