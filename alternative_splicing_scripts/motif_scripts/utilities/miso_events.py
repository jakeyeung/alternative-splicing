'''
Created on 2014-01-17

@author: jyeung
'''

import sys
import os
import csv


def get_inclusion_exclusion(incl_file=None, excl_file=None):
    '''
    Read inclusion and exclusion fasta file containing miso events,
    return dictionary containing 
    
    {miso_event: 'inclusion'|'exclusion' ,...}
    '''
    # Define keys
    incl_str = 'inclusion'
    excl_str = 'exclusion'
    outdic = {}
    
    for jfile, incl_or_excl in zip([incl_file, excl_file], 
                                   [incl_str, excl_str]):
        with open(jfile, 'rb') as readfile:
            for line in readfile:
                if line.startswith('>'):    # contains miso event
                    # remove first character '>', strip \n
                    miso_event = line[1:].strip()
                    if miso_event not in outdic:
                        outdic[miso_event] = incl_or_excl
                    else:
                        print 'Duplicate in miso event: %s' %miso_event
                        print 'Exiting...'
                        sys.exit()
    return outdic

def go_up_two_dirs(tomtom_path):
    '''
    Given tomtom path, which should be two directories deep from
    the directory containing region information,
    extract the region.
    '''
    return os.path.basename(os.path.dirname(os.path.dirname(tomtom_path)))

def get_tomtom_hits(meme_dir, rel_path):
    '''
    Given meme dir, retrieve all the tomtom hits, returning a dictionary
    of the form:
    
    Other inputs:
    rel_path: Relative path (from meme_dir) to matched RBPs 
    
    gene1, gene2 matching to motif 1 inclusion
    {Motif 1 inclusion: [gene1, gene2]}
    '''
    # Define column names in RBP file
    motif_colname = '#Query ID'    # motif number
    rbp_colname = 'Target ID'    # CSV containing gene name
    
    # change directory to meme_dir, may solve long path issues
    os.chdir(meme_dir)
    
    # init outdic
    outdic = {}
    
    # create paths to tomtom results
    alldirs = os.listdir(os.getcwd())
    meme_dirs = []    # dirs containing meme results
    for jdir in alldirs:
        if jdir.startswith('exon') or jdir.startswith('intron'):
            meme_dirs.append(jdir)
    # append rel_path to meme_dir to create absolute paths to files
    meme_paths = []
    for jdir in meme_dirs:
        path = os.path.join(os.getcwd(), jdir, rel_path)
        if os.path.exists(path):
            meme_paths.append(path)
    # Iterate each path, store info into outdic 
    for path in meme_paths:
        with open(path, 'rb') as readfile:
            jreader = csv.reader(readfile, delimiter='\t')
            try:
                header = jreader.next()
            except StopIteration:    # empty file
                continue
            for row in jreader:
                motif_numb = row[header.index(motif_colname)]
                motif_name = row[header.index(rbp_colname)]    # CSV
                # first element in CSV is gene
                gene = motif_name.split(',')[0]
                # create key for dic combining position, inclusion, motif numb
                region = go_up_two_dirs(path)
                # replace underscore with spaces for region
                region = ' '.join(region.split('_'))
                # append motif number to region
                # Key should match plot_meme_motif_summary.py
                key = ' '.join(['Motif', motif_numb, region])
                # init subdic with empty list if it doesnt exist
                if key not in outdic:
                    outdic[key] = []
                if gene not in outdic[key]:    # no redundant genes
                    outdic[key].append(gene)
    return outdic
                