'''
Created on 2014-01-01

@author: jyeung

This script is to be run after
summarize_meme_results.py, which creates and stores
an output dic.

It reads the pickle dic from summarize_meme_results.py
and returns GERP conservation scores associated with
its genome coordinates. 

The general strategy is to index chromo:position with
RS scores. Then read the dic to get RS scores across
a range of start-ends, returning average RS scores across
each range.
'''

import sys
import os
import csv
from optparse import OptionParser
import pickle
from multiprocessing import Process, Queue
from utilities.gerp_utilities import get_chr_list
from summarize_meme_results import \
    get_dic_subkeys, write_outdic_to_file

def init_gerp_dic(motif_dic):
    '''
    Read pickle containing motif dic, then 
    extract all the chromosomes and genomic coords
    that will need to be read in RS score.
    
    motif dic is expected to be of form:
    {miso_event: {region: {subkey1: subval1 ...}}}
    
    Returns:
    gerp_dic of form:
    {chr1: {coord1: None, coord1: None...}
    chrY: {...}} 
    '''
    # Define subkey constants. Should match get_dic_subkeys()
    genomic_coord_str = 'genomic_coordinate'
    
    # Init chromosome list
    chr_list = get_chr_list()

    # Init gerp_dic
    gerp_dic = {}
    
    # Create empty subdics, by chromosmoe
    for chromosome in chr_list:
        gerp_dic[chromosome] = {}
    
    # Iterate through the different subdics, retrieving genomic coords
    # store genomic coords to gerp_dic
    coord_count = 0
    for miso_event in motif_dic:
        for region in motif_dic[miso_event]:
            # get colon separated chr_start_end value.
            chr_start_end_list = \
                motif_dic[miso_event][region][genomic_coord_str]
            # loop through list of coords (mostly length 1)
            for chr_start_end in chr_start_end_list:
                if chr_start_end is not None:
                    # Split by colon
                    chr_start_end_split = chr_start_end.split(':')
                    # Get chromo, start, ends separately.
                    # Start is inclusive in RS score file.
                    # End is not inclusive (so we subtract by 1).
                    chromo = chr_start_end_split[0]
                    start = int(chr_start_end_split[1])
                    end = int(chr_start_end_split[2])
                    # Update gerp_dic containing empty dics
                    # with coordinates as keys, from start to end
                    # Note: end is not inclusive.
                    for coordinate in range(start, end):
                        gerp_dic[chromo][coordinate] = None
                        coord_count += 1
    print '%s coordinates stored into gerp_dic.' %coord_count
    return gerp_dic

def get_gerp_fpath(chromosome, gerp_dir):
    '''
    Given chromosome and gerp directory, return full path
    to open the RS scores file.
    '''
    rs_scores_fname = '.'.join([chromosome, 'maf', 'rates'])
    rs_scores_fpath = os.path.join(gerp_dir, rs_scores_fname)
    return rs_scores_fpath

def add_rs_scores_to_gerp_dic(gerp_dic, chromo, gerp_dir, queue_obj):
    '''
    ***THIS IS USED WITH MULTIPROCESSING*** That's why queue_obj.put()
    is there instead of returning a dic.
    
    Given gerp_dic of form:
    gerp_dic of form:
    {chr1: {coord1: None, coord1: None...}
    chrY: {...}} 
    
    open gerp_file from gerp_dir for input chromo, 
    get rs_scores for relevant positions (by looking at subkeys),
    update subvalues for each position with rs_scores.
    
    Put gerp_dic[chromo] into queue_obj. This is important 
    for multiprocessing because
    when we update, we want to update only gerp_dic[chromo].
    
    Updating gerp_dic will result in "Nones" in a lot of RS scores
    that are not in respective chromosome.  
    '''
    
    # Get gerp fpath
    gerp_fpath = get_gerp_fpath(chromo, gerp_dir)
    
    '''
    Find max in coordinates, so we know when to
    stop iterating through rows when we've foudn all the
    RS scores we need.
    
    Also, if gerp_dic[chromo] is an empty dic, just return
    gerp_dic. This chromosome has no relevant coordinates
    '''
    coords_list = gerp_dic[chromo].keys()
    if len(coords_list) == 0:
        print 'Empty coordinates list in chromosome: %s\nSkipping...\n' %chromo
        queue_obj.put((gerp_dic[chromo], chromo))
        return
    else:
        max_coord = max(coords_list)
    
    '''
    open file, iterate through lines. Keep track
    of lines we iterate. The row number is equivalent to 
    genomic coordinate.
    '''
    update_count = 0
    with open(gerp_fpath, 'rb') as gerp_file:
        rs_scores_reader = csv.reader(gerp_file, delimiter='\t')
        for coordinate, row in enumerate(rs_scores_reader):
            if coordinate in gerp_dic[chromo]:
                # coordinate matches coordinate in gerpdic, 
                # retrieve its RS score.
                rs_score = float(row[1])
                # Update coordinat gerp dic with rs score
                gerp_dic[chromo][coordinate] = rs_score
                update_count += 1
            elif coordinate > max_coord:
                break    # no need to iterate, we've got all our RS scores.
            else:
                continue    # didnt match chromo, but should still iterate.
    print '%s scores extracted for chromosome: %s' \
        %(update_count, chromo)  
    queue_obj.put((gerp_dic[chromo], chromo))    # for multiprocessing
    return

def get_rs_score_subkey():
    rs_score_subkey = 'avg_rs_score'
    return rs_score_subkey

def get_chromo_start_end(genomic_coordinate):
    '''
    Given string, colon separated chromo:start:end,
    return chromo, start, end.
    '''
    coord_split = genomic_coordinate.split(':')
    chromo = coord_split[0]
    start = int(coord_split[1])
    end = int(coord_split[2])
    return chromo, start, end 

def get_avg_rs_score(genomic_coordinate, gerp_dic):
    '''
    Given genomic coordinate (chromo:start:end),
    find gerp values from gerp_dic.
    '''
    chromo, start, end = get_chromo_start_end(genomic_coordinate)
    rs_scores_list = []
    for coord in range(start, end):
        rs_scores_list.append(gerp_dic[chromo][coord])
    avg_rs_score = float(sum(rs_scores_list)) / len(rs_scores_list)
    return avg_rs_score

def update_motif_dic_with_gerp_dic(motif_dic, gerp_dic):
    '''
    Updates motif dic with RS scores found in GERP dic
    
    I expect sub-sub-key called genomic_coordinate to exist.
    '''
    coord_subkey = 'genomic_coordinate'
    avg_score_subkey = get_rs_score_subkey()
    
    for miso_event in motif_dic:
        for region in motif_dic[miso_event]:
            try:
                coords_list = \
                    motif_dic[miso_event][region][coord_subkey]
            except KeyError:
                print 'Could not find subkey: %s in %s' \
                    %(coord_subkey, motif_dic[miso_event][region])
            # Calculate average RS score from gerp dic given chromo:start:end
            avg_rs_scores_list = []
            for coord in coords_list:
                avg_rs_scores_list.append(get_avg_rs_score(coord, gerp_dic))
            # add avg_rs_score to motif_dic
            if avg_score_subkey not in motif_dic[miso_event][region]: 
                motif_dic[miso_event][region][avg_score_subkey] = \
                    avg_rs_scores_list
            else:
                print 'Error: did not expect duplicate in %s' \
                    %motif_dic[miso_event][region]
    return motif_dic

def main():
    usage = 'usage: %prog pickle_filepath gerp_directory output_file\n'\
        'Two args must be specified in commandline: \n'\
        '1) Path to pickle from summarize_meme_results.py\n'\
        '2) Directory containing GERP RS score text files by chromosome\n'\
        '3) Output file to which results will be written.\n'
    parser = OptionParser(usage=usage)
    parser.add_option('-p', '--gerp_pickle_fname', dest='gerp_pickle_fname',
                      default='gerp_pickle.pkl',
                      help='gerp scores pickle filename.\n'\
                        'Default gerp_pickle.pkl')
    parser.add_option('-l', '--presaved_gerp_dic_path', dest='gerp_presaved_path',
                      default=None,
                      help='If a gerp dic has been presaved, use this flag to'\
                        'indicate the file path to directly open the gerp file. '\
                        'Reduces need for multiprocessing.')
    parser.add_option('-o', '--updated_dic_fname', dest='updated_dic_fname',
                      default='meme_summary.gerp_updated.pkl',
                      help='Filename to updated meme_summary.pkl'
                        ' with gerp scores,\n'\
                        'Default "meme_summary.gerp_updated.pkl"')
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print 'Incorrect number of parameters specified.'
        print usage
        sys.exit()
    motif_pickle_path = args[0]
    gerp_dir = args[1]
    output_path = args[2]
    
    # parse options
    gerp_pickle_fname = options.gerp_pickle_fname
    gerp_presaved_pkl_path = options.gerp_presaved_path
    updated_dic_fname = options.updated_dic_fname
    
    # Load motif dic, obtained from summarize_meme_results
    pickle_file = open(motif_pickle_path, 'rb')
    motif_dic = pickle.load(pickle_file)
    pickle_file.close()
    
    # Get chromosome list
    chr_list = get_chr_list()
    
    # Create dic containing chromosomes and genomic coordinates
    # relevant to our motifs by reading pickled dictionary.
    gerp_dic = init_gerp_dic(motif_dic)
    
    # Multithread only if presaved_gerp_dic_path flag is not None.
    if gerp_presaved_pkl_path == None:
        # BEGIN: MULTITHREADING
        print 'Beginning multiprocessing.'
        q = Queue()
        process_list = []
        # For each chromosome, open relevant gerp file and retrieve
        # RS scores associated with the coordinates.
        for chromosome in chr_list:
            print 'Calculating RS scores for %s' %chromosome
            p = Process(target=add_rs_scores_to_gerp_dic,
                        args=(gerp_dic, chromosome, gerp_dir, q))
            process_list.append(p)
            p.start()
        
        for _ in chr_list:
            # Update dic for every process that started.
            # the actual chromosoem doesn't matter. It's the
            # number of iterations that matter.
            (gerp_dic_chromo, chromo) = q.get()
            # Find which chromosome this came from by looking at key
            gerp_dic[chromo].update(gerp_dic_chromo)
            
        # Wait for all threads to be done before continuing.
        for p in process_list:
            p.join()
        
        # END: MULTITHREADING
        print 'Done multiprocessing.'
        
        # save gerp dic as pickle
        gerp_pickle_dir = os.path.dirname(motif_pickle_path)
        gerp_pickle_fpath = os.path.join(gerp_pickle_dir, gerp_pickle_fname)
        with open(gerp_pickle_fpath, 'wb') as gerp_pickle_file:
            pickle.dump(gerp_dic, gerp_pickle_file, -1)
        print 'Saved pickle object to: %s' %gerp_pickle_fpath 
    
    else:
        # Try to open pickle path from options.gerp_path
        with open(gerp_presaved_pkl_path, 'rb') as presaved_pkl:
            pickle.load(presaved_pkl)
        print 'Loaded presaved pickle from %s' %gerp_presaved_pkl_path
            
    # Update motif_dic with gerp_dic
    motif_dic = update_motif_dic_with_gerp_dic(motif_dic, gerp_dic)
    
    # Save updated motif dic to pickle
    updated_dic_dir = os.path.dirname(motif_pickle_path)
    updated_dic_path = os.path.join(updated_dic_dir, updated_dic_fname)
    with open(updated_dic_path, 'wb') as updated_dic_file:
        pickle.dump(motif_dic, updated_dic_file, -1)
    print 'Updated dic saved to: %s' %updated_dic_path
    
    # Write updated gerp_dic to file
    # add GERP scores as a subkey in subkey_list
    subkeys_list = get_dic_subkeys()
    rs_score_subkey = get_rs_score_subkey()
    subkeys_list.append(rs_score_subkey)
    
    # Write updated motif dic to file
    write_outdic_to_file(motif_dic, output_path, subkeys_list)

if __name__ == '__main__':
    main()