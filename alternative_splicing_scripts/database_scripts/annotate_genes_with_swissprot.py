'''
Created on 2013-12-18

@author: jyeung

Requires biopython!
'''

import sys
import csv
from Bio import SwissProt
from optparse import OptionParser
from create_dna_protein_summary_file import get_subkeys

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
    
def main():
    usage = 'usage: %prog [opt] protein_summary_file output_filename'\
        '\nTwo arguments must be specified in command line:\n'\
        '1) Protein summary file (output from create_dna_protein_summary_file.py)\n'\
        '2) Output file'
    parser = OptionParser(usage=usage)
    
    # colnames for lfq data
    parser.add_option('-f', '--swiss_prot_flat_file', dest='swissprot_file',
                      default='G:/jyeung/projects/alternative_splicing/input/uniprot_database/uniprot_sprot.dat',
                      help='File path to swiss prot flat file')
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print 'Not enough args specified.\n%s' %usage
        sys.exit()
    
    candidate_genes_file = args[0]
    outfile = args[1]
    swissprot_file = options.swissprot_file
    
    # Index candidate genes for matching to swissprot.
    genes_dic, rowcount = index_candidate_genes(candidate_genes_file)
    print '%s amino acid sequences indexed.' %rowcount
    
    # define genename identifier
    name_str = 'Name'
    
    for count, record in enumerate(SwissProt.parse(open(swissprot_file))):
        print record.organism
        print record.gene_name.split(';')
        print record.features
        print record.sequence
        raw_input()
    print 'itereated through %s swissprot records.' %count
    
if __name__ == '__main__':
    main()