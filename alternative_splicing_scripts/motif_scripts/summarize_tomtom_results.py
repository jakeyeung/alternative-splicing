'''
Created on 2014-01-08

@author: jyeung

Read every tomtom file in directory, summarize it into one textfile.

Written with cassette exons in mind.
'''

import sys
import os
import csv
from optparse import OptionParser
from summarize_meme_results import \
    str_to_boolean, get_region_of_interest_from_filepath
    
def get_sites_from_line(line):
    '''
    Given
    MOTIF  2    width =   15   sites =  20   llr = 306   E-value = 3.9e-027
    Return 20
    '''
    # Double split to get n_sites
    n_sites = line.split("sites =")[1].split('   ')[0]
    return int(n_sites)

def create_motif_key(motif_number):
    '''
    Given motif number, create motif key (used for dictionaries)
    motif_number an integer.
    Returns motif key
    '''
    motif_key = '_'.join(['motif', str(motif_number)])
    return motif_key

def get_n_sites_from_meme(meme_file):
    '''
    Read meme text file, return dic of form:
    {motif_1: {n_sites: 27} ...}
    '''
    # define colnames
    n_sites_str, _, _, _, _ = get_subkeys()
    
    # init dic
    n_sites_dic = {}
    
    motif_count = 1    # first motif is motif_1
    # build motif search string, double space join.
    searchstring = '  '.join(['MOTIF', str(motif_count)])
    with open(meme_file, 'rb') as jfile:
        for line in jfile:
            if line.startswith(searchstring):
                # get sites from line
                n_sites = get_sites_from_line(line)
                # create keyname: e.g. motif_1
                motif_key = create_motif_key(motif_count)
                if motif_key not in n_sites_dic:
                    n_sites_dic[motif_key] = {n_sites_str: n_sites}
                else:
                    print 'Error: duplicate motif_key for %s' %motif_key
                    sys.exit()
                motif_count += 1
                searchstring = '  '.join(['MOTIF', str(motif_count)])
    return n_sites_dic

def get_tomtom_results(tomtom_file):
    '''
    Read TomTom file, create dic of form:
    {motif_1: {gene: geneA, p-value: X, E-value: Y, q-value: Z}}
    #Query ID    Target ID    Optimal offset    p-value    E-value    q-value
    '''
    # define colnames
    pval_colname = 'p-value'
    eval_colname = 'E-value'
    qval_colname = 'q-value'
    motif_colname = '#Query ID'
    gene_colname = 'Target ID'
    
    # define subkey strings
    # prevent confusion by saying gene, not Target ID for subkey
    _, gene_str, pval_str, eval_str, qval_str = get_subkeys()
    
    # init outdic
    tomtom_dic = {}
    
    with open(tomtom_file, 'rb') as myfile:
        myreader = csv.reader(myfile, delimiter='\t')
        try:
            header = myreader.next()
        except StopIteration:    # empty file
            return {}
        for row in myreader:
            motif = row[header.index(motif_colname)]
            gene = row[header.index(gene_colname)]
            p_val = row[header.index(pval_colname)]
            e_val = row[header.index(eval_colname)]
            q_val = row[header.index(qval_colname)]
            # Create key, write to dic
            motif_key = create_motif_key(motif)
            # If motif_key not in tomtom, init subkey with empty lists
            if motif_key not in tomtom_dic:
                tomtom_dic[motif_key] = {}
                for subkey in [pval_str, eval_str, 
                               qval_str, gene_str]:
                    tomtom_dic[motif_key][subkey] = []
            # Append subvals to list in subkey
            for subkey, subval in \
                zip([pval_str, eval_str, qval_str, gene_str], 
                    [p_val, e_val, q_val, gene]):
                # try to save as float, otherwise just append without change.
                try:
                    tomtom_dic[motif_key][subkey].append(float(subval))
                except ValueError:
                    tomtom_dic[motif_key][subkey].append(subval)
    return tomtom_dic

def get_subkeys():
    '''
    Get column names used for output or for dictionary keys.
    '''
    n_sites_str = 'n_sites'
    gene_str = 'gene'
    pval_str = 'pval'
    eval_str = 'eval'
    qval_str = 'qval'
    return n_sites_str, gene_str, pval_str, eval_str, qval_str
    
def main():
    usage = 'usage: %prog meme_results_file output_file\n'\
        'Two args must be specified in commandline: \n'\
        '1) Directory containing meme results (at exon/intron level).\n'\
        '2) Output file to which results will be written.\n'
    
    parser = OptionParser(usage=usage)
    
    parser.add_option('-f', '--tomtom_filename', dest='tomtom_filename',
                      default='rbp_matches/candidate_rbps.txt',
                      help='tomtom output filename (built from arg 1). '\
                        'Default rbp_matches/candidate_rbps.txt')
    parser.add_option('-m', '--meme_filename', dest='meme_filename',
                      default='meme.txt',
                      help='Meme output filename. Default meme.txt')
    parser.add_option('-p', '--pickle_filename', dest='pickle_filename',
                      default='tomtom_summary.pkl',
                      help='Output filename for pickled dictionary.'\
                      'Default tomtom_summary.pkl')
    parser.add_option('-n', '--null_mode', dest='null_mode',
                      default=False,
                      help='Option flag to retrieve entire fasta genomic \n'\
                        'coordinate rather than motif genomic coordinate.\n'\
                        'Useful for obtaining a null distribution of gerp '\
                        'scores.')
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print 'Incorrect number of parameters specified.'
        print usage
        sys.exit()
    meme_html_dir = args[0]
    output_path = args[1]
    
    # parse tomtom/meme_filename options
    tomtom_filename = options.tomtom_filename
    meme_filename = options.meme_filename
    # parse null mode options
    null_mode = options.null_mode
    null_mode = str_to_boolean(null_mode)
    
    # Define column names
    n_sites_str, gene_str, pval_str, eval_str, qval_str = get_subkeys()
    motif_str = 'motif_id'
    region_str = 'region'
    subkeys = [gene_str, region_str, motif_str, n_sites_str, pval_str, 
               eval_str, qval_str]
    
    # create list of tomtom paths linking to tomtom files
    # for each intronic and exonic region.
    all_dirs = os.listdir(meme_html_dir)
    paths_tup_list = []
    
    for mydir in all_dirs:
        if mydir.startswith('intron_') or mydir.startswith('exon_'):
            tomtom_path = os.path.join(meme_html_dir, mydir, tomtom_filename)
            meme_path = os.path.join(meme_html_dir, mydir, meme_filename)
            # only append if it is a true path to a file
            if os.path.isfile(tomtom_path) and os.path.isfile(meme_path):
                # Save as tuple, meme first, tomtom second.
                paths_tup_list.append((meme_path, tomtom_path))
    
    # Itereate through each region, get number of sites per motif and 
    # candidate RBPs for each motif. Save it all to a dic.
    outdic = {}
    for meme_file, tomtom_file in paths_tup_list:
        region = get_region_of_interest_from_filepath(meme_file, 
                                                      incl_excl=True)
        outdic[region] = {}
        # create dic of {motif_1: {n_sites: 10}}
        n_sites_dic = get_n_sites_from_meme(meme_file)
        outdic[region].update(n_sites_dic)
        # get tomtom results
        tomtom_dic = get_tomtom_results(tomtom_file)
        # update dic to include tomtom info
        # use try in case tomtom doesn't have info on that motif.
        for motif in outdic[region]:
            try:
                outdic[region][motif].update(tomtom_dic[motif])
            except KeyError:
                pass
    # Add motif information to outdic by concatenating region with motif
    for region in outdic:
        for motif in outdic[region]:
            # create motif_id from region and motif
            motif_id = ':'.join([region, motif])
            outdic[region][motif][motif_str] = motif_id
    
    # Write to output, ready for ggplot2.
    writecount = 0
    with open(output_path, 'wb') as writefile:
        jwriter = csv.writer(writefile, delimiter='\t')
        # write header, include motif id to end of subkeys. 
        jwriter.writerow(subkeys)
        for region in outdic:
            for motif in outdic[region]:
                # Only include subkeys containing gene names
                if gene_str not in outdic[region][motif]:
                    continue
                # one row per gene hit (for ggplot2)
                subdic = outdic[region][motif]    # reduces typing later
                # itereate lists in parallel
                for gene, pval, e_val, qval in zip(subdic[gene_str], 
                                                  subdic[pval_str], 
                                                  subdic[eval_str], 
                                                  subdic[qval_str]):
                    # order must match subkeys() here...
                    writerow = []
                    for val in [gene, region, subdic[motif_str], subdic[n_sites_str], 
                                pval, e_val, qval]:
                        writerow.append(val)
                    jwriter.writerow(writerow)
                    writecount += 1
    print '%s rows written to: %s' %(writecount, output_path)
                
if __name__ == '__main__':
    main()