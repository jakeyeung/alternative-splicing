'''
Created on 2013-12-14

@author: jyeung

Given a discovered motif in an intron region, can we say
if those events that make up the motif also have 
other motifs in other intron regions?
'''

import sys
import os
from optparse import OptionParser

def get_seq_names_from_meme_results(meme_path, motif_numb):
    '''
    Read meme path, knowing motif number, find the sequence names
    taht match to that motif.
    
    Output list of seq names
    
    Must match string twice. First time it matches for log likelihood
    matrix. Second time it matches for sequence names.
    '''
    # define a match string to know when you've reached the
    # motif you want...
    matchstring = 'BL   MOTIF %s' %motif_numb
    # end_of_seq_names = '//'    # begins row where there are no more seq names.
    
    seq_names_list = []
    matched_to_seq_names_row = False
    with open(meme_path, 'rb') as readfile:
        for myline in readfile:
            if myline.startswith(matchstring):
                if matched_to_seq_names_row == False:
                    # Not yet at seq names, but next
                    # time it matches it wil be seq names.
                    # So switch to True then continue.

                    matched_to_seq_names_row = True
                    continue
                else:
                    '''
                    # We're now at row containnig seq names!
                    # Go to next row then
                    # let's start collecting our seq names
                    We stop when the lines start with //, meaning it is end.
                    '''
                    seq_name_full_row = readfile.next()
                    while not seq_name_full_row.startswith('//'):
                        # Split seq name by double space, grab first element.
                        seq_name = seq_name_full_row.split('  ')[0]
                        seq_names_list.append(seq_name)
                        seq_name_full_row = readfile.next()
                    # loop is done, get out. No need to iterate rows.
                    break
    return seq_names_list
    
def main():
    usage = 'usage: %prog meme_results motif_number '\
        'meme_results_to_compare output_file'\
        '\n'\
        'Four args must be specified in commandline: \n'\
        '1) Meme HTML file containing meme results.\n'\
        '2) Integer indicating motif number to match.'\
        ' Must not exceed number of motifs in meme result.\n'\
        '3) MEME html from from another region in which to find matches.\n'\
        '-h to display this help message.\n'
    parser = OptionParser(usage=usage)
    (_, args) = parser.parse_args()
    if len(args) < 3:
        print 'Incorrect number of parameters specified.'
        print usage
        sys.exit()
    meme_html_path = args[0]
    motif_number = args[1]
    meme_comparison_path = args[2]
    # output_path = args[3]
    
    # Read meme results, get sequence names from specified motif.
    seq_names_list = get_seq_names_from_meme_results(meme_html_path, 
                                                    motif_number)
    n_seq_names = len(seq_names_list)
    intron_region = os.path.basename(os.path.dirname(meme_html_path))
    print '%s sequence names matched to %s motif %s' %(n_seq_names,
                                                       intron_region, 
                                                       motif_number)
    
    # Read meme comparison results, get sequences from all motifs, count
    # up from 1 until no more motifs are matched 
    # (motif_numb_contains_seqs False)
    motif_number_contains_sequences = True    # init
    motif_compare_number = 1    # initialize
    while motif_number_contains_sequences:
        seq_names_compare_list = \
            get_seq_names_from_meme_results(meme_comparison_path, 
                                            motif_compare_number)
        if len(seq_names_compare_list) == 0:
            motif_number_contains_sequences = False
        else:
            match_count = 0
            for compare_seq in seq_names_compare_list:
                if compare_seq in seq_names_list:
                    match_count += 1
            meme_dir = os.path.basename(os.path.dirname(meme_comparison_path))
            print 'Motif %s from %s match count: %s/%s' %(motif_compare_number, 
                                                       meme_dir, 
                                                       match_count,
                                                       n_seq_names)
            motif_compare_number += 1
if __name__ == '__main__':
    main()