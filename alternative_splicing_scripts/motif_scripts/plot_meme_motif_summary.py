'''
Created on 2014-01-03

@author: jyeung

requires matplotlib

Plots conservation of motifs. 

Requires running of three scripts: 

1) summarize_meme_results.py 
2) summarize_meme_results_with_gerp_scores.py
3) append_gene_names_to_textfile.py

After running that, we will try to determine the
distribution of conservation for these motifs.
'''

import sys
import csv
import os
from optparse import OptionParser
import matplotlib.pyplot as plt
from utilities import miso_events, plot_functions

def get_regions(meme_header):
    '''
    In meme header, column names may be colon separated.
    Example:
    
    intron_1_3p:motif_relative_end
    intron_2_3p:motif_sequence
    
    Parse meme_header, find colon separated colnames, and
    return the list of regions in the header.
    
    We do this because we want this script to run on
    MXE as well.
    
    '''
    regions_list = []
    for colname in meme_header:
        colname_split = colname.split(':')
        if len(colname_split) != 2:
            # if length != 2, means this colname is not
            # of form region:feature
            continue
        else:
            # Get first element from split, which should
            # be the region, e.g. intron_1_3p in
            # motif_relative_end
            regions_list.append(colname_split[0])
    return list(set(regions_list))

def get_avg_rs_scores_as_dic():
    pass
    
def create_tomtom_key(motif_id, region):
    '''
    Given region and motif id, recreate a key
    that will allow access to tomtom
    
    motif_id looks like: 
    Motif 2 inclusion
    
    region looks like:
    intron_2_3p
    
    We want to return a tomtom key...
    tomtom keys look lke:
    Motif 2 intron 2 5p inclusion
    '''
    # split into lists
    region = region.split('_')
    motif_id = motif_id.split(' ')
    
    # insert region into second 3rd position in list (index 2)
    # we iterate the region so that we don't get a list within a list.
    for piece in region:
        motif_id.insert(-1, piece)
    
    # Join with spaces
    return ' '.join(motif_id)

def main():
    usage = 'usage: %prog meme_gerp_genename_filepath output_filepath\n'\
        'Requires two input arguments:\n'\
        '1) textfile output from '\
            'summarize_meme_results_with_gerp_scores\n'\
        '2) inclusion fasta file\n'\
        '3) exclusion fasta file\n'\
        '4) meme dir containing meme results'
    parser = OptionParser(usage=usage)    
    (_, args) = parser.parse_args()
    
    if len(args) < 5:
        print 'Four arguments need to be specified in command line.\n'
        print usage
        sys.exit()
    meme_summarypath = args[0]
    incl_fasta = args[1]
    excl_fasta = args[2]
    meme_dir = args[3]
    
    # define column name suffix (string after colon in colname)
    gerp_str = 'avg_rs_score'
    motif_numb_str = 'motif_number'
    # miso event has no col name suffix, this is entire colname
    miso_colname = 'miso_event'
    # define rel path to tomtom files from meme dir
    rel_path = os.path.join('rbp_matches', 'candidate_rbps.txt')
    # define plot title
    mytitle = 'GERP Score Comparison: Hits vs Non-Hits'
    
    # get dictionary containing inclusion and exclusion for miso event
    incl_excl_dic = miso_events.get_inclusion_exclusion(incl_file=incl_fasta, 
                                                        excl_file=excl_fasta)
    tomtom_dic = miso_events.get_tomtom_hits(meme_dir, rel_path)
    region_gerp_scores = {}    # gerp scores, indexed by region.
    with open(meme_summarypath, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        header = myreader.next()
        regions = get_regions(header)
        # init output dic with empty lists
        for region in regions:
            region_gerp_scores[region] = {}
        
        for row in myreader:
            # get gerp score in each region.
            # beware of empty values.
            for region in regions:
                subdic = region_gerp_scores[region]
                gerp_colname = ':'.join([region, gerp_str])
                motif_numb_colname = ':'.join([region, motif_numb_str])
                gerp_score = row[header.index(gerp_colname)]
                motif_numb = row[header.index(motif_numb_colname)]
                miso_event = row[header.index(miso_colname)]
                incl_or_excl = incl_excl_dic[miso_event]
                motif_id = ' '.join(['Motif', motif_numb, incl_or_excl])
                if gerp_score is not '':
                    if motif_id not in subdic:
                        subdic[motif_id] = []
                    subdic[motif_id].append(float(gerp_score))          
                    
    avg_scores_in_tomtom = []
    avg_scores_not_in_tomtom = []
    # Plot histogram of average scores.
    for region in region_gerp_scores:
        for motif_id in region_gerp_scores[region]:
            tomtom_key = create_tomtom_key(motif_id, region)
            if tomtom_key in tomtom_dic:
                avg_scores_in_tomtom += region_gerp_scores[region][motif_id]
            else:
                avg_scores_not_in_tomtom += region_gerp_scores[region][motif_id]
    
    conserved_counts_in_tomtom = 0
    conserved_counts_not_in_tomtom = 0
    for s in avg_scores_in_tomtom:
        if s >= 2:
            conserved_counts_in_tomtom += 1
    for s in avg_scores_not_in_tomtom:
        if s >= 2:
            conserved_counts_not_in_tomtom += 1
    
    for avg_scores, mylabel in zip([avg_scores_in_tomtom, 
                                    avg_scores_not_in_tomtom], 
                                   ['Motif with matching RBPs', 
                                    'Motif without matching RBPs']):         
        plot_functions.plot_density(avg_scores, 
                                    mytitle=mytitle, 
                                    mylabel=mylabel)
    plt.legend()
    plt.show()

if __name__ == '__main__':
    main()