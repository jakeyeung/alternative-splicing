'''
Created on 2013-12-22

@author: jyeung

Parse meme html output and retrieve 
start and stops of each motif hit.

Currently only works for miso events from cassette exons!

'''

import sys
import os
import csv
import random
from optparse import OptionParser
import pickle

def get_motif_start_end(motif_line):
    '''
    Motif line example:
    chr22:19964938:19965109:-@chr22:19964229:19964246:-@chr22:19963209:19963280:-  ( 23) GGGAGGGCATGGGGG 1
    
    Split by spaces, retrieve the 3rd element from split'd line, 
    which would be "23)", then remove the character from string, ).
    
    Return 23, which should be an integer.
    '''
    motif_linesplit = motif_line.split(' ')
    motif = motif_linesplit[4]
    try:
        motif_start = int(motif_linesplit[3][:-1])
    except ValueError:
        print 'Expected %s to be an integer.' %motif_start
        raw_input()
    motif_end = motif_start + len(motif)
    try:
        return int(motif_start), int(motif_end)
    except ValueError:
        print 'Could not extract motif start from: %s' %motif_line
        print '%s must be an integer.' %motif_start
        sys.exit()
        
def get_sequence_from_motif_line(motif_line):
    '''
    Motif line example:
    chr22:19964938:19965109:-@chr22:19964229:19964246:-@chr22:19963209:19963280:-  ( 23) GGGAGGGCATGGGGG 1
    
    Split by spaces, retrieve the 5th element from split line.
    Return GGGAGGGCATGGGGG     
    '''
    return motif_line.split(' ')[4]
        
def get_motif_name_from_motif_line(motif_line):
    '''
    Motif line example:
    chr22:19964938:19965109:-@chr22:19964229:19964246:-@chr22:19963209:19963280:-  ( 23) GGGAGGGCATGGGGG 1
    
    Split by spaces, retrieve the 1st element from split line.
    Return chr22:19964938:19965109:-@chr22:19964229:19964246:-@chr22:19963209:19963280:-
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
    
    Example 3:
    chr8:144679518:144679845:-@chr8:144676045:144676241:-@chr8:144674885:144675112:-
    if intron_number == 1, intron_3p_or_5p == "5p", seq_length = 100
    Return:
    start = 144679418 (144679518 - 100)
    end = 144679518
    
    Intron number should be int
    '''
    strand = get_strand_from_miso_event(miso_event)
    intron_starts, intron_ends = get_intron_starts_ends_list(miso_event)
    
    # get seq length - 1 because start coordinates are inclusive.
    seq_length_adj = seq_length - 1
    
    # get intron start and ends corresponding to intron index
    intron_index = intron_number - 1

    if strand == '+':
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
    elif strand == '-':
        # Get seq length at 3p or 5p site of intron start/end
        if intron_3p_or_5p == '5p':
            intron_end = intron_ends[intron_index]
            intron_start = intron_end - seq_length_adj
        elif intron_3p_or_5p == '3p':
            intron_start = intron_starts[intron_index]
            intron_end = intron_start + seq_length_adj
        else:
            print 'Expected intron_3p_or_5p to be "5p" or "3p". %s found.' \
                %intron_3p_or_5p
            sys.exit()
    return intron_start, intron_end
        
def get_exon_starts_ends(miso_event, exon_number, seq_length):
    '''
    Given:
    chr22:19964938:19965109:-@chr22:19964229:19964246:-@chr22:19963209:19963280:-
    
    if exon_number == 2:
    Return 
    exon_start = 19964229
    exon_end = 19964246
    '''
    exon_index = exon_number - 1    # starts at 0
    exon_coords = miso_event.split('@')[exon_index].split(':')
    exon_start = int(exon_coords[1])
    exon_end = int(exon_coords[2])
    return exon_start, exon_end
    
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
        exon_number = int(region_split[1])
        start, end = get_exon_starts_ends(miso_event, exon_number, strand)
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
    else:
        print 'Expected %s to begin with "exon" or "intron". %s found.' \
            %(region_of_interest, exon_or_intron)
        sys.exit()
    return chromo, start, end

def get_region_of_interest_from_filepath(filepath, incl_excl=False):
    '''
    I expect filepath to be for example:
    ~.../fasta_unshuffled/intron_1_3p_inclusion/meme.html
    
    We want intron_1_3p to be returned if incl_excl == False
    We want intron_1_3p_inclusion if incl_excl == True
    '''
    # Find gene region (intron or exon) based on the dir name of html file
    myregion = os.path.basename(os.path.dirname(filepath))
    
    # myregion still has _inclusion, let's remove that.
    myregion_split = myregion.split('_')
    # check first element in split is exon or intron, throw error if not.
    if myregion_split[0] != 'exon' and myregion_split[0] != 'intron':
        print 'Expected exon or intron in string: %s' %myregion_split[0]
        sys.exit()
    # rejoin, removing the last element (incl/excl) in list if 
    # incl_excl == False. Do not remove anything if incl_excl == True.
    if incl_excl == False:
        return '_'.join(myregion_split[:-1])
    elif incl_excl == True:
        return '_'.join(myregion_split)
    else:
        print 'incl_excl must be True or False. %s found.' %incl_excl

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
                        print 'Overwriting...'
                    relevant_line = readfile.next()
    return seq_lengths_dic

def get_motif_start_end_coordinates(seq_start, seq_end, 
                                    motif_rel_start, motif_rel_end, 
                                    strand):
    '''
    Purpose:
    Given sequence start and ends and motif relative starts/ends, return
    motif starts and stops relative to genomic coordinates.
    Must take strand into consideration, because sequence start/end is
    based on UCSC coordinates (i.e. start > end, always, even if negative 
    strand)
    Inputs:
    seq_start: genomic start coordinate (chromosome agnostic)
    seq_end: genomic end coordinate
    motif_rel_start: when the motif begins, relative to input sequence
        (motif_rel_start = 0 is the first basepair in sequence, also
        means motif_rel_start == seq_start if motif_rel_start == 0)
    motif_rel_end: when the motif ends, relative to input sequence
    strand: "+" or "-".
    
    Outputs:
    motif_start_coord: genomic start coordinate of motif
    motif_end_coord: genomic end coordinate of motif
    
    motif_start_coord > motif_end_coord always, no matter what strand.
    This helps visualization in UCSC.
    '''
    # motif_rel_end must be deincremented by 1 because coords are inclusive
    motif_rel_end_adj = motif_rel_end - 1
    if strand == '+':
        motif_start_coord = seq_start + motif_rel_start
        motif_end_coord = seq_start + motif_rel_end_adj
    elif strand == '-':
        '''
        If strand is -, seq_end is actually beginning of input sequence
        because seq_start > seq_end always, as a rule.
        
        So even though seq_end - motif_rel_start is begining of motif
        relative to transcription, we name it as end coord because it is
        the farther down if we take + strand as reference.
        '''
        motif_end_coord = seq_end - motif_rel_start
        motif_start_coord = seq_end - motif_rel_end_adj
    else:
        print 'Expected %s to be "+" or "-"' %strand
        sys.exit()
    return motif_start_coord, motif_end_coord

def init_dic(mydic, keys_list):
    '''
    Given a dictionary, create a key:value where
    key is in keys_list and value is empty list.
    
    Note: only create an empty list IF key does not exist in dic
    '''
    for key in keys_list:
        if key not in mydic:
            mydic[key] = []
    return mydic

def get_csv_from_list(mylist):
    '''
    Given a list, return csv. 
    All elements in list gets converted to string.
    '''
    mylist = [str(i) for i in mylist]
    return ','.join(mylist)

def get_tomtom_motifs(tomtom_dic, meme_html_path):
    '''
    Inputs:
        tomtom_dic: contains motif matches given region
        meme_html_path: can be parsed to get region
        
    Output:
        list of motif numbers that correspond to a tomtom match.
    '''
    region = \
        get_region_of_interest_from_filepath(meme_html_path, incl_excl=True)
    motif_matches = tomtom_dic[region]['motif_number']
    if len(motif_matches) > 1:
        motif_matches.sort()
    # order in ascending order, might be important.
    return motif_matches

def store_info_to_output_dic(outdic, miso_event, region_of_interest, 
                             subkeys_list, motif_info_list):
    '''
    Do two checks: one if miso_event is initialized,
    second if region_of_interest
    Rationale:
        if we loop across many meme outputs, we will
        need dictionary to hold a summary of miso events
        across all regions of interest.
        
    ALSO do one more check:
    if in TomTom mode, make sure current motif number
    is in list of motif numbers matching to TomTom.
    If it does not, then do not write to dic.
    '''
    if miso_event not in outdic:
        # intiialize if miso_event not yet a key
        outdic[miso_event] = {}
        outdic[miso_event][region_of_interest] = {}
    elif region_of_interest not in outdic[miso_event]:
        # initialize only region of interest if not yet subkey
        outdic[miso_event][region_of_interest] = {}
    else:
        pass
    # create key:empty lists
    outdic[miso_event][region_of_interest] = \
        init_dic(outdic[miso_event][region_of_interest], 
                 subkeys_list)
    '''
    # Add motif information to output dic in lists
    # in case multiple motifs matching to same event
    '''
    for subkey, subval in zip(subkeys_list, motif_info_list):
        # append subval to subkey list
        outdic[miso_event][region_of_interest][subkey].\
            append(subval)

def update_dic_with_motifs(outdic, meme_html_path, region_of_interest, 
                           seq_lengths_dic, 
                           null_mode=False,
                           null_length=9,
                           tomtom_mode=False,
                           tomtom_dic=None,
                           use_entire_intron=False):
    '''
    Loop through meme html file, for every motif in meme html file, retrieve
    all the miso events that match the motif as well as
    the genomic coordinates of the motif and relative motif position with 
    respect to the region.
    
    null_mode: option flag to not retrieve motif start and ends, but rather
    only the fasta seq start and ends. This is useful if you want to 
    do some null distribution of gerp scores.
    
    Null mode: two things change: relative start/end and motif start/end.
    Relative start/end is randomly chosen intronic region
    rel_start: random number btwn 0 to seq_length-null_length
    rel_end: rel_start + null_length
    
    Tomtom mode: if tomtom mode is True, tomtom_dic cannot be None.
    We will check motif to see if motif matches a tomtom motif, if it
    does not, then we skip it.
    '''
    # Get subkeys list
    subkeys_list = get_dic_subkeys()
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
    
    # if tomtom mode, get a list of motifs that correspond to tomtom motifs
    if tomtom_mode:
        if tomtom_dic is None:
            print 'Error: In Tomtom mode: Tomtom dic cannot be None.'
            sys.exit()
        tomtom_motif_numbers = get_tomtom_motifs(tomtom_dic, meme_html_path)
    
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
                    # BEGIN: RETRIEVE MOTIF INFORMATION
                    miso_event = get_motif_name_from_motif_line(motif_line)
                    strand = get_strand_from_miso_event(miso_event)
                    # motif start/end, relative to beginning of fasta sequence
                    if not null_mode:
                        motif_rel_start, motif_rel_end = \
                            get_motif_start_end(motif_line)
                    else:
                        if use_entire_intron:
                            motif_rel_start = 0
                            motif_rel_end = seq_lengths_dic[miso_event]
                        else:
                            # In null mode: randomly pick number as
                            # start (bounded by seq_length).
                            # end is start + null_length
                            lower_bnd = 0
                            #TODO: upper_bnd cannot be negative!
                            upper_bnd = seq_lengths_dic[miso_event] - null_length + 1
                            motif_rel_start = \
                                random.randrange(lower_bnd, upper_bnd)
                            motif_rel_end = motif_rel_start + null_length
                            
                    chromo, seq_start, seq_end = \
                        get_seq_start_end_from_miso_event(miso_event, 
                                                          region_of_interest,
                                                          seq_lengths_dic)
                    # Check seq_end > seq_start, or else error.
                    if seq_start > seq_end:
                        print 'Error: expected end(%s) > start(%s).'
                        sys.exit()
                    # Get motif start and end genomic coordinates.
                    motif_start, motif_end = \
                        get_motif_start_end_coordinates(seq_start, 
                                                        seq_end,
                                                        motif_rel_start,
                                                        motif_rel_end, 
                                                        strand)
                    
                    # Concatenate to get genomic coordinate with chromo
                    genomic_coord = ':'.join([chromo, str(motif_start), 
                                              str(motif_end)])
                    # Get motif sequence
                    motif_seq = get_sequence_from_motif_line(motif_line)
                    # END: RETRIEVE MOTIF INFORMATION
                    
                    # BEGIN: Create list of relevant information ready 
                    # to be stored in dic.
                    '''
                    subkeys_list is in order:
                        1) motif_rel_start
                        2) motif_rel_end
                        3) genomic_coordinate
                        4) motif_sequence
                        5) motif_number
                        6) Avg RS score (conservation)
                    '''
                    motif_info_list = [motif_rel_start, motif_rel_end, 
                                       genomic_coord, motif_seq, 
                                       motif_number]
                    # BEGIN: store information to output dic
                    if tomtom_mode:
                        if motif_number in tomtom_motif_numbers:
                            store_info_to_output_dic(outdic, miso_event, 
                                                     region_of_interest, 
                                                     subkeys_list, 
                                                     motif_info_list)
                    else:
                        store_info_to_output_dic(outdic, miso_event, 
                                                 region_of_interest, 
                                                 subkeys_list, 
                                                 motif_info_list)
                    # END: store information to output dic
                    motif_line = readfile.next()
                # Out of while loop, means we should increment motif number 
                motif_number += 1
                search_motif_string = ''.join([searchstring, 
                                               str(motif_number)])
    return outdic

def get_dic_subkeys():
    '''
    Get predefined dictionary subkeys.
    '''
    # Define subkeys used within outdic
    motif_rel_start_str = 'motif_relative_start'
    motif_rel_end_str = 'motif_relative_end'
    genomic_coord_str = 'genomic_coordinate'
    sequence_str = 'motif_sequence'
    motif_number_str = 'motif_number'
    subkeys_list = [motif_rel_start_str, motif_rel_end_str, 
                    genomic_coord_str, sequence_str, 
                    motif_number_str]
    return subkeys_list

def get_gerp_filepath():
    '''
    Get gerp filepath
    '''
    gerp_filepath = \
        'G:\jyeung\projects\alternative_splicing\input'\
        '\gerp_conservation_score\gerp_scores_by_chr'
    return gerp_filepath

def save_dic_as_pickle(mydic, picklepath):
    '''
    Given a dictionary and path to file, 
    check if file exists. If does not exist,
    save dic as pickle to that path.
    If it does exist, return warning message.
    '''
    if os.path.isfile(picklepath):
        print 'Pickle path: %s already exists.'\
        'Overwriting...' %picklepath
    pickle_output = open(picklepath, 'wb')
    pickle.dump(mydic, pickle_output, -1)    # highest protocol
    pickle_output.close()
    print 'Pickle object saved to: %s' %picklepath

def write_outdic_to_file(outdic, output_path, subkeys_list):
    '''
    Now that we have our wonderful dictionary. let's write it to a 
    textfile. This needs to be general enough so that
    if we add anything to our outdic, we will still be able to
    use this function to write to file.
    '''
    
    # Get regions list from outdic:
    regions_list = []
    for miso_event in outdic:
        for region in outdic[miso_event]:
            if region not in regions_list:
                regions_list.append(region)
    regions_list = sorted(regions_list)
    
    # Create output column names, append region (eg intron_1_3p to subkey)
    miso_event_str = 'miso_event'
    # Initialize column name list
    output_colnames = [miso_event_str]
    # append region to subkey
    for region in regions_list:
        for subkey in subkeys_list:
            output_colnames.append(':'.join([region, subkey]))
    
    # Write to output file
    with open(output_path, 'wb') as outfile:
        outwriter = csv.writer(outfile, delimiter='\t')
        # write colnames
        outwriter.writerow(output_colnames)
        for writecount, miso_event in enumerate(outdic):
            # init empty list and add miso event to list
            row_to_write = [miso_event]
            for region in regions_list:
                # subkey might not exist in dic, skip if it doesn't exist
                if region not in outdic[miso_event]:
                    for subkey in subkeys_list:
                        row_to_write.append(None)
                else:
                    for subkey in subkeys_list:
                        # Get value_list from subkey
                        subval_list = outdic[miso_event][region][subkey]
                        # Collapse list to CSV, append to row_to_write
                        subval_csv = get_csv_from_list(subval_list)
                        row_to_write.append(subval_csv)
            outwriter.writerow(row_to_write)
    print '%s rows written to %s:' %(writecount, output_path)
    
def str_to_boolean(mystr):
    '''
    Parse string, return True or False as boolean.
    '''
    if mystr in ['True', 'true', 'TRUE', True, 'T', 't']:
        return True
    elif mystr in ['False', 'false', 'FALSE', False, 'F', 'f']:
        return False
    else:
        print 'Error: expected %s to be True or False.' %mystr
        return None
    
def get_motif_numb_from_id(motif_id):
    '''
    Input a string like: intron_2_3p_inclusion:motif_1
    Return 1, as an integer.
    '''
    motif_str = motif_id.split(':')[1]
    motif_numb = motif_str.split('_')[1]
    return int(motif_numb)
    
def get_tomtom_dic(tomtom_path):
    '''
    Read tomtom summary, retrieve dictionary of form:
    {intron_2_3p_inclusion: {region: intron_2_3p, 
                            incl_excl = inclusion, 
                            motif_number: [1]},...,}
    '''
    tomtom_dic = {}
    with open(tomtom_path, 'rb') as readfile:
        jreader = csv.reader(readfile, delimiter='\t')
        header = jreader.next()
        for row in jreader:
            region = row[header.index('region')]    # does not hv incl/excl
            incl_excl = row[header.index('incl_or_excl')]
            # combine region with incl_excl to get region_incl
            full_region = '_'.join([region, incl_excl])
            # get motif from motif_id
            motif_id = row[header.index('motif_id')]
            motif_numb = get_motif_numb_from_id(motif_id)
            # update tomtom if not yet updated
            if full_region not in tomtom_dic:
                tomtom_dic[full_region] = {}
                tomtom_dic[full_region]['region'] = region
                tomtom_dic[full_region]['incl_excl'] = incl_excl
                tomtom_dic[full_region]['motif_number'] = []
            # Append motif number to list if not in list
            if motif_numb not in tomtom_dic[full_region]['motif_number']:
                tomtom_dic[full_region]['motif_number'].append(motif_numb)
    return tomtom_dic
            
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
    parser.add_option('-e', '--include_exons', dest='include_exons',
                      default=False,
                      help='True or False: option flag to include exon motifs '\
                      'in analysis or not.')
    parser.add_option('-f', '--meme_filename', dest='meme_filename',
                      default='meme.html',
                      help='Meme output filename. Default meme.html')
    parser.add_option('-p', '--pickle_filename', dest='pickle_filename',
                      default='meme_summary.pkl',
                      help='Output filename for pickled dictionary.'\
                      'Default meme_summary.pkl')
    parser.add_option('-n', '--null_mode', dest='null_mode',
                      default=False,
                      help='Option flag to retrieve entire fasta genomic \n'\
                        'coordinate rather than motif genomic coordinate.\n'\
                        'Useful for obtaining a null distribution of gerp '\
                        'scores.')
    parser.add_option('-t', '--tomtom_summarypath', dest='tomtom_path',
                      default=None,
                      help='Reads from tomtom summary (summarize_tomtom_results.py).\n'\
                        'Summarizes meme results only for hits in TomTom summary.\n'\
                        'Input a filepath to tomtom summary to activate this option.')
    parser.add_option('--null_length', dest='null_length',
                      default=9,
                      help='Length of "null" motifs. Default is 9.')
    parser.add_option('-A', '--use_entire_intron', dest='use_entire_intron',
                      action="store_true", default=False)
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print 'Incorrect number of parameters specified.'
        print usage
        sys.exit()
    meme_html_dir = args[0]
    output_path = args[1]
    
    # parse include_exons option
    include_exons = options.include_exons
    include_exons = str_to_boolean(include_exons)
    # parse meme_filename options
    meme_filename = options.meme_filename
    # parse pickle filename option
    pickle_filename = options.pickle_filename
    # parse null mode options
    null_mode = options.null_mode
    null_mode = str_to_boolean(null_mode)
    # parse tomtom option
    tomtom_path = options.tomtom_path
    if tomtom_path is None:
        tomtom_mode = False
    else:
        tomtom_mode = True
    # parse null length option
    null_length = int(options.null_length)
    
    # Tell user if they are in null mode
    if null_mode:
        print 'Processing in NULL mode '\
            '(takes a random region of length %s)' %null_length
        # Set random seed for getting random intronic region
        random.seed(0)
        
    # Tell user if they are in tomtom mode, then get tomtom dic
    if tomtom_mode:
        print 'Only considering motifs that match tomtom summary.'
        print 'Tomtom summary file: %s' %tomtom_path
        tomtom_dic = get_tomtom_dic(tomtom_path)
    else:
        tomtom_dic = None
    
    # create list of meme html paths linking to meme.html files
    # for each intronic and exonic region.
    all_dirs = os.listdir(meme_html_dir)
    meme_paths = []
    # include only dirs starting wih intron_ and maybe exon_ 
    # depending on option flags
    for mydir in all_dirs:
        # If in TomTom mode, chec, if dir is in tomtom dic, if not then
        # ignore it. This limits paths to only ones matching to tomtom motifs
        if tomtom_mode:
            if mydir not in tomtom_dic:
                continue
        if include_exons == False:
            if mydir.startswith('intron_'):
                # append path to dir
                mypath = os.path.join(meme_html_dir, mydir, meme_filename)
                # only append if it is a true path to a file
                if os.path.isfile(mypath):
                    meme_paths.append(mypath)
        elif include_exons == True:
            if mydir.startswith('intron_') or mydir.startswith('exon_'):
                mypath = os.path.join(meme_html_dir, mydir, meme_filename)
                # only append if it is a true path to a file
                if os.path.isfile(mypath):
                    meme_paths.append(mypath)
        else:
            print 'include_exons must be True or False. %s found.' \
                %include_exons
            sys.exit()
    
    # Loop through meme html paths
    regions_list = []
    # Create output dic to store information
    outdic = {}
    for meme_html_path in meme_paths:
        print 'Retrieving motif info from %s' %meme_html_path
        
        # get region of interest from filename
        region_of_interest = \
            get_region_of_interest_from_filepath(meme_html_path)
        # store region of interest for later output
        regions_list.append(region_of_interest)
        
        # create dictionary of miso_event: sequence length from meme html file
        seq_lengths_dic = get_seq_lengths_from_meme_file(meme_html_path)
        
        # Fill dic with motif info
        outdic = update_dic_with_motifs(outdic, meme_html_path, 
                                        region_of_interest, seq_lengths_dic,
                                        null_mode=null_mode,
                                        null_length=null_length,
                                        tomtom_mode=tomtom_mode,
                                        tomtom_dic=tomtom_dic,
                                        use_entire_intron=options.use_entire_intron)
    
    subkeys_list = get_dic_subkeys()
    write_outdic_to_file(outdic, output_path, subkeys_list)
    
    # Save outdic to pickle file. In same directory as output file
    outdir = os.path.dirname(output_path)
    picklepath = os.path.join(outdir, pickle_filename)
    save_dic_as_pickle(outdic, picklepath)
    
if __name__ == '__main__':
    main()