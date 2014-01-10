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
    n_sites_str, _, _, _, _, _ = get_subkeys()
    
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

def get_incl_or_excl(region):
    '''
    From region, determine if it is inclusion or exclusion.
    Example of region:
    intron_1_3p_inclusion,
    Return inclusion in that case.
    '''
    return region.split('_')[-1]

def get_tomtom_results(tomtom_file, region):
    '''
    Read TomTom file, create dic of form:
    {motif_1: {gene: geneA, p-value: X, E-value: Y, 
        q-value: Z, motifid: A, incl_or_excl: incl/excl}}
    #Query ID    Target ID    Optimal offset    p-value    E-value    q-value
    
    We need region in input arg to determine if it is inclusion ore xclusion.
    '''
    # define colnames
    pval_colname = 'p-value'
    eval_colname = 'E-value'
    qval_colname = 'q-value'
    motif_colname = '#Query ID'
    gene_colname = 'Target ID'
    
    # define subkey strings
    # prevent confusion by saying gene, not Target ID for subkey
    _, gene_str, pval_str, eval_str, qval_str, incl_or_excl_str = \
        get_subkeys()
    
    # init outdic
    tomtom_dic = {}
    gene_dic = {}
    
    # Determine if it is inclusion or exclusion from region
    incl_or_excl = get_incl_or_excl(region)
    
    with open(tomtom_file, 'rb') as myfile:
        myreader = csv.reader(myfile, delimiter='\t')
        try:
            header = myreader.next()
        except StopIteration:    # empty file
            return {}
        for row in myreader:
            motif = row[header.index(motif_colname)]
            # convert ELAVL3,M232_0.6,I -> ELAVL3
            gene = row[header.index(gene_colname)].split(',')[0]
            p_val = row[header.index(pval_colname)]
            e_val = row[header.index(eval_colname)]
            q_val = row[header.index(qval_colname)]
            # Create key, write to dic
            motif_key = create_motif_key(motif)
            motif_gene_key = ':'.join([motif_key, gene])
            if motif_gene_key not in gene_dic:
                gene_dic[motif_gene_key] = {}
                for subkey in [pval_str, eval_str,
                               qval_str, gene_str,
                               incl_or_excl_str]:
                    gene_dic[motif_gene_key][subkey] = []
            # Append subvals to list in subkey in gene dic
            for subkey, subval in \
                zip([pval_str, eval_str, qval_str, gene_str, incl_or_excl_str], 
                    [p_val, e_val, q_val, gene, incl_or_excl]):
                # try to save as float, otherwise just append without change.
                try:
                    gene_dic[motif_gene_key][subkey].append(float(subval))
                except ValueError:
                    gene_dic[motif_gene_key][subkey].append(subval)
    # Add gene_dic values to tomtom_dic, but take average
    # significance values in string
    for motif_gene_key in gene_dic:
        motif = motif_gene_key.split(':')[0]
        gene = motif_gene_key.split(':')[1]
        # If motif_key not in tomtom, init subkey with empty lists
        if motif not in tomtom_dic:
            tomtom_dic[motif] = {}
            for subkey in [pval_str, eval_str, 
                           qval_str, gene_str, 
                           incl_or_excl_str]:
                tomtom_dic[motif][subkey] = []
        for subkey in [pval_str, eval_str, qval_str]:
            val_list = gene_dic[motif_gene_key][subkey]
            avg_val = sum(val_list) / float(len(val_list))
            tomtom_dic[motif][subkey].append(avg_val)
        # Add incl or exclusion and gene name to tomtom_dic
        # Append gene name to tomtom dic
        gene_list = list(set(gene_dic[motif_gene_key][gene_str]))
        incl_excl_list = list(set(gene_dic[motif_gene_key][incl_or_excl_str]))
        for subkey, mylist in zip([gene_str, incl_or_excl_str], 
                                  [gene_list, incl_excl_list]):
            # check length of set is 1, append
            if len(mylist) == 1:
                tomtom_dic[motif][subkey].append(mylist[0])
            else:
                print 'expected list length to be 1: %s' %mylist
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
    incl_or_excl_str = 'incl_or_excl'
    return n_sites_str, gene_str, pval_str, \
            eval_str, qval_str, incl_or_excl_str
    
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
    n_sites_str, gene_str, pval_str, eval_str, qval_str, incl_or_excl_str = \
        get_subkeys()
    motif_str = 'motif_id'
    region_str = 'region'
    subkeys = [gene_str, region_str, motif_str, n_sites_str, pval_str, 
               eval_str, qval_str, incl_or_excl_str]
    
    # create list of tomtom paths linking to tomtom files
    # for each intronic and exonic region.
    all_dirs = os.listdir(meme_html_dir)
    paths_tup_list = []
    
    for mydir in all_dirs:
        if mydir.startswith('intron_') or mydir.startswith('exon_'):
            tomtom_path = os.path.join(meme_html_dir, mydir, tomtom_filename)
            meme_path = os.path.join(meme_html_dir, mydir, meme_filename)
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
        tomtom_dic = get_tomtom_results(tomtom_file, region)
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
                for gene, pval, e_val, qval, incl_excl in zip(subdic[gene_str], 
                                                              subdic[pval_str], 
                                                              subdic[eval_str], 
                                                              subdic[qval_str],
                                                              subdic[incl_or_excl_str]):
                    # order must match subkeys() here...
                    # set region to remove inclusion/exclusion
                    region_no_inclexcl = '_'.join(region.split('_')[:-1])
                    writerow = []
                    for val in [gene, region_no_inclexcl, subdic[motif_str], subdic[n_sites_str], 
                                pval, e_val, qval, incl_excl]:
                        writerow.append(val)
                    jwriter.writerow(writerow)
                    writecount += 1
    print '%s rows written to: %s' %(writecount, output_path)
                
if __name__ == '__main__':
    main()