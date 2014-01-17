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
from optparse import OptionParser
import matplotlib.pyplot as plt
from utilities import miso_events

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

def plot_histogram(values_list, n_bins, mytitle, mylabel):
    '''
    Given list of values, plot histogram.
    '''
    plt.hist(values_list, n_bins, histtype='step', 
             stacked=False, fill=True, alpha=0.5,
             label=mylabel)
    plt.title(mytitle)

def main():
    usage = 'usage: %prog meme_gerp_genename_filepath output_filepath\n'\
        'Requires two input arguments:\n'\
        '1) textfile output from '\
            'summarize_meme_results_with_gerp_scores\n'\
        '2) inclusion fasta file\n'\
        '3) exclusion fasta file\n'
    parser = OptionParser(usage=usage)    
    (_, args) = parser.parse_args()
    
    if len(args) < 3:
        print 'Three arguments need to be specified in command line.\n'
        print usage
        sys.exit()
    meme_summarypath = args[0]
    incl_fasta = args[1]
    excl_fasta = args[2]
    
    # define column name suffix (string after colon in colname)
    gerp_str = 'avg_rs_score'
    motif_numb_str = 'motif_number'
    # miso event has no col name suffix, this is entire colname
    miso_colname = 'miso_event'
    
    # get dictionary containing inclusion and exclusion for miso event
    incl_excl_dic = miso_events.get_inclusion_exclusion(incl_file=incl_fasta, 
                                                        excl_file=excl_fasta)
    print incl_excl_dic
    
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

    # Plot histogram of average scores.
    for region in region_gerp_scores:
        for motif_id in region_gerp_scores[region]:
            avg_scores = region_gerp_scores[region][motif_id]
            plot_histogram(avg_scores, n_bins=25, 
                           mytitle=region, mylabel=motif_id)
        plt.legend()
        plt.show()

if __name__ == '__main__':
    main()