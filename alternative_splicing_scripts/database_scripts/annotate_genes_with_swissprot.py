'''
Created on 2013-12-18

@author: jyeung

Requires biopython!
'''

import sys
import csv
import pickle
from Bio import SwissProt
from optparse import OptionParser
from create_dna_protein_summary_file import get_subkeys
from index_unitprot_db import get_uniprot_subkeys

def index_candidate_genes(gene_file):
    '''
    Read gene file (created from create_dna_protein_summary_file.py)
    Returns a dictionary of indexed gene information.
    {gene_name: {amino_acid_seq:[], miso_event:[], reading_frame:[]}}
    '''
    gene_dic = {}
    
    # get same subkeys as the subkeys used to create summary file
    reading_frame_colname, \
    gene_name_colname, \
    _, \
    _, \
    _, \
    _ = get_subkeys()
    
    # get other relevant column names
    miso_event_colname = 'miso_event'
    amino_acid_seq_colname = 'amino_acid_sequence'
    
    # index gene file to dic
    with open(gene_file, 'rb') as readfile:
        genereader = csv.reader(readfile, delimiter='\t')
        geneheader = genereader.next()
        for rowcount, row in enumerate(genereader):
            # get relevant row information
            gene_name = row[geneheader.index(gene_name_colname)]
            amino_acid_seq = row[geneheader.index(amino_acid_seq_colname)]
            miso_event = row[geneheader.index(miso_event_colname)]
            reading_frame = row[geneheader.index(reading_frame_colname)]
            
            '''
            # store row info to dic, gene name as key
            # initialize empty subdic and empty lists if
            # gene name key is not yet initialized...
            '''
            if gene_name not in gene_dic:
                gene_dic[gene_name] = {}
                for subkey in [amino_acid_seq_colname, 
                               miso_event_colname, 
                               reading_frame_colname]:
                    gene_dic[gene_name][subkey] = []
            else:    # gene name already exists in dic, so do nothing
                pass
            '''
            With properly initialized dics, we now have either subkeys with 
            empty lists or subkeys with non-emtpy lists, we will append
            new data into those lists.
            '''
            for subkey, subvalue in zip([amino_acid_seq_colname, 
                                         miso_event_colname, 
                                         reading_frame_colname], 
                                        [amino_acid_seq, 
                                         miso_event, 
                                         reading_frame]):
                gene_dic[gene_name][subkey].append(subvalue)
    return gene_dic, rowcount

def extract_gene_name_from_record(record):
    '''
    From swissprot record, extract gene name...
    '''
    try:
        gene_name = record.gene_name.split(';')[0].split('=')[1]
    except IndexError:
        print 'Could not find gene name in %s' %record.gene_name
        raw_input()
        gene_name = None
    return gene_name

def get_amino_acid_start_end(exon_aa_seq, gene_aa_seq):
    '''
    Given amino acid of entire gene, search where the exon amino acid
    lies within the gene.
    Return start and stop of amino acid alignment. 
    I expect 100% alignment, so no fancy algorithms here, just straight
    searching.
    '''
    if exon_aa_seq in gene_aa_seq:
        # find when exon_aa starts within gene_aa
        # index retrieves start, add length for end.
        exon_aa_start = gene_aa_seq.index(exon_aa_seq)
        exon_aa_end = exon_aa_start + len(exon_aa_seq)
        return exon_aa_start, exon_aa_end
    else:
        return None, None
    
def append_dic_if_feature_within_start_end(exon_start, exon_end, 
                                           amino_acid_seq,
                                           uniprot_dic, gene_key, feature,
                                           output_dic):
    '''
    Given start and end, check if a particular feature within
    a gene inside uniprot_dic matches start and ends in the
    feature annotation.
    
    Return all instances where it matches in a dictionary object.
    
    amino acid sequence comes from a particular exon.
    
    Dictionary format:
    {amino_acid_sequence: {feature: {[start], [end], [description]}}}
    '''
    # get uniprot subkeys for accessing feature starts, stops and descriptions
    start_subkey, end_subkey, descript_subkey = get_uniprot_subkeys()
    # define two additional subkeys: exon_start and exon_end
    exon_start_subkey = 'exon_start'
    exon_end_subkey = 'exon_end'
    
    # initialize match_count
    match_count = 0
    
    # get start, end, description from uniprot dic
    feature_start_list = uniprot_dic[gene_key][feature][start_subkey]
    feature_end_list = uniprot_dic[gene_key][feature][end_subkey]
    descript_list = uniprot_dic[gene_key][feature][descript_subkey]

    '''
    # iterate feature start/end in parallel, ask if it is within 
    # the exon start/end range.
    Criteria for if it is NOT within range is:
    exon_start > feature_end
    exon_end < feature_start
    '''
    for feature_start, feature_end, descript in zip(feature_start_list, 
                                                    feature_end_list,
                                                    descript_list):
        if exon_start > feature_end or exon_end < feature_start:
            # feature outside of relevant range, go to next start/end
            continue
        else:
            # feature within relevant range, store to output dic
            # intialize relevant keynames if not yet initialized already.
            output_keyname = amino_acid_seq
            if output_keyname not in output_dic:
                output_dic[output_keyname] = {}
                output_dic[output_keyname][feature] = {}
                for subkey in [start_subkey, end_subkey, descript_subkey, 
                               exon_start_subkey, exon_end_subkey]:
                    output_dic[output_keyname][feature][subkey] = []
            else:
                # If amino acid seq already exists in output_dic, just return
                # output dic (i.e. go to next amino acid seq)
                return output_dic, match_count
            # store values into subkey
            for subkey, subval in \
                zip([start_subkey, end_subkey, descript_subkey, 
                     exon_start_subkey, exon_end_subkey],
                    [feature_start, feature_end, descript, 
                     exon_start, exon_end]):
                output_dic[output_keyname][feature][subkey].append(subval)
            match_count += 1
    return output_dic, match_count

def annotate_exon_seq(exon_aa_seq, gene_aa_seq, uniprot_dic, gene_key, 
                      output_dic):
    '''
    Given exon seq, gene seq, uniprot dic and gene key, 
    Update output_dic to include uniprot annotations.
    '''
    # define amino acid seq subkey, important because it's not a feature.
    amino_acid_seq_subkey = 'amino_acid_sequence'
    
    # define match count
    match_count = 0
    match_count_total = 0
    
    # Begin annotating exon seq for uniprot features...            
    aa_start, aa_end = get_amino_acid_start_end(exon_aa_seq, 
                                                gene_aa_seq)
    # if searching failed, start and end would return None
    if aa_start != None and aa_end != None:
        '''
        # knowing start and ends, search uniprot_dic for modified residues
        # iterate through subkeys (except amino_acid_seq, since it's 
        # not a feature), check if that feature is relevant within
        # our exon.
        '''
        for feature in uniprot_dic[gene_key]:
            if feature == amino_acid_seq_subkey:
                # this is not a feature, but an AA sequence.
                # so skip to next feature.
                continue
            output_dic, match_count = \
                append_dic_if_feature_within_start_end(aa_start, 
                                                       aa_end, 
                                                       exon_aa_seq,
                                                       uniprot_dic, 
                                                       gene_key, 
                                                       feature, 
                                                       output_dic)
            match_count_total += match_count
    return output_dic, match_count_total
    
def main():
    usage = 'usage: %prog [opt] protein_summary_file output_filename'\
        '\nTwo arguments must be specified in command line:\n'\
        '1) Protein summary file '\
        '(output from create_dna_protein_summary_file.py)\n'\
        '2) Output file'
    parser = OptionParser(usage=usage)
    
    # colnames for lfq data
    parser.add_option('-f', '--uniprot_pickle_dic', dest='uniprot_pkl_path',
                      default='G:/jyeung/projects/alternative_splicing/'\
                      'input/uniprot_database/uniprot_sprot.with_seq.pkl',
                      help='File path to swiss prot flat file')
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print 'Not enough args specified.\n%s' %usage
        sys.exit()
    
    candidate_genes_file = args[0]
    outfile = args[1]
    uniprot_pkl_path = options.uniprot_pkl_path
    
    # Index candidate genes for matching to swissprot.
    genes_dic, rowcount = index_candidate_genes(candidate_genes_file)
    print '%s amino acid sequences indexed.' %rowcount
    
    # Load pickled dictionary to file
    uniprot_dic = pickle.load(open(uniprot_pkl_path, 'rb'))
    print 'Pickled dic loaded from %s' %uniprot_pkl_path
    
    # define genename identifier
    name_str = 'Name'
    
    # define genes_dic subkeys
    amino_acid_seq_subkey = 'amino_acid_sequence'
    
    # initialize an output dic, which will store matched exon annotations
    output_dic = {}
    
    # init counter
    total_matches = 0
    
    '''
    # begin iterating through candidate genes and checking through uniprot
    # database if they have relevant features WITHIN the amino_acid seq of 
    # interest
    '''
    for gene_count, gene in enumerate(genes_dic):
        '''
        Iterate through genes, assemble gene_key that
        matches uniprot_dic, find incl/excl exons's amino acid
        amino acid number (through some type of searching),
        then grab all annotations surrounding incl/excl sequence
        from uniprot.
        '''
        # 1) Assemble gene key, e.g. geneA -> Name=geneA
        gene_key = '='.join([name_str, gene])
        # if gene_key not found in uniprot_dic, go to next gene
        if gene_key not in uniprot_dic:
            continue
        # 2) Get amino acid number by searching exon sequence
        # with gene sequence in uniprot.
        # gene_aa_seq is a string, exon_aa_seqs is list of strings.
        gene_aa_seq = uniprot_dic[gene_key][amino_acid_seq_subkey]    # string
        exon_aa_seqs = genes_dic[gene][amino_acid_seq_subkey]    # is a list
        # Get non-redundant list
        exon_aa_seqs = list(set(exon_aa_seqs))
        # iterate over non-redundant list to find matches
        for exon_aa_seq in exon_aa_seqs:
            output_dic, match_count = \
                annotate_exon_seq(exon_aa_seq, 
                                  gene_aa_seq, 
                                  uniprot_dic, 
                                  gene_key,
                                  output_dic)
            total_matches += match_count
    print '%s matches after searching for %s genes.' %(total_matches, 
                                                       gene_count)
    
    # Write protein annotations to output file by recreating
    # summary file and appending uniprot annotations to the end.
    
    for key in output_dic:
        print key
        print output_dic[key]
        raw_input()
    
if __name__ == '__main__':
    main()
