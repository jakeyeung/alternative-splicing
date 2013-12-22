'''
Created on 2013-12-19

@author: jyeung
UnitProt contains info from across many organisms.
Create a pickle dictionary, reducing the dat file to only
homo sapiens and relevant information of the protein (e.g. modified
residues)
'''

import sys
import os
from optparse import OptionParser
import pickle
from Bio import SwissProt

def get_uniprot_subkeys():
    '''
    Returns subkeys used to identify start, end, description
    of features.
    '''
    # def subkeys
    start_str = 'start'
    end_str = 'end'
    description_str = 'description'
    return start_str, end_str, description_str

def index_uniprot_database(uniprot_file, species='homosapiens'):
    '''
    Index uniprot database for homosapien entries
    Use Name=genename as key (don't want to use genename as key
    because duplicate keys don't work in python, and I may have
    a dictionary with genename as key)
    '''
    
    # def constants
    if species == 'homosapiens':
        organism_str = 'Homo sapiens (Human).'
    else:
        print 'Other organisms not yet coded in this script...'
        print 'Only homosapiens works at the moment.'
        sys.exit()
    
    # def subkeys
    start_subkey, end_subkey, description_subkey = get_uniprot_subkeys()
    amino_acid_seq_subkey = 'amino_acid_sequence'
        
    uniprot_dic = {}
    for count, record in enumerate(SwissProt.parse(open(uniprot_file))):
        if record.organism != organism_str:
            pass
        else:    # organism matches, we want this in our dic.
            '''
            # get key, Name=genename
            # but record.gene_name could have four possibilities:
            1) Synonyms
            2) Ordered locus names
            3) ORF names
            4) Name
            We want only Name, try to iterate through list, and grab
            only the entry containing Name=
            '''
            gene_name_list = record.gene_name.split(';')
            # iterate the list, and find teh one that says Name=
            for gene_name in gene_name_list:
                if gene_name.startswith('Name='):
                    gene_name_key = gene_name
                    break    # we found it, break out of for loop
            '''
            # Check if gene name is alerady in dic, if it is not,
            # then initialize empty lists otherwise, append to 
            '''
            if gene_name_key not in uniprot_dic:
                # not in dic, so let's init
                uniprot_dic[gene_name_key] = {}
            '''
            # we do not know what the subkeys are actually, we will take them 
            # all from record.features, which is a list of tuples with 
            (keyname, from, to, description)
            '''
            aa_sequence = record.sequence
            # Add sequence to subdic
            uniprot_dic[gene_name_key][amino_acid_seq_subkey] = aa_sequence
            for feature in record.features:
                keyname = feature[0]
                residue_start = feature[1]
                residue_end = feature[2]
                description = feature[3]
                
                if keyname not in uniprot_dic[gene_name_key]:
                    # initialize keyname with dic, then subkeys with list
                    uniprot_dic[gene_name_key][keyname] = {}
                    
                    # if subkey not exist, init empty list in a subkey
                    for subkey in start_subkey, end_subkey, description_subkey:
                        uniprot_dic[gene_name_key][keyname][subkey] = []
                else:    # list exists, we will append to it.
                    pass
                
                for subkey, subvalue in zip([start_subkey, 
                                             end_subkey, 
                                             description_subkey], 
                                            [residue_start, 
                                             residue_end, 
                                             description]):
                    uniprot_dic[gene_name_key][keyname][subkey].append(subvalue)
    return uniprot_dic, count

def main():
    usage = 'usage: %prog [opt] protein_summary_file output_filename'\
        '\nOne argument must be specified in command line:\n'\
        '1) Output of dictionary pickle'
    parser = OptionParser(usage=usage)
    
    parser.add_option('-f', '--uniprot_flat_file', dest='uniprot_file',
                      default='G:/jyeung/projects/alternative_splicing/input/uniprot_database/uniprot_sprot.dat',
                      help='File path to swiss prot flat file')
    (options, args) = parser.parse_args()
    if len(args) < 1:
        print 'Not enough args specified.\n%s' %usage
        sys.exit()
    
    output_pickle_path = args[0]
    uniprot_file = options.uniprot_file
    
    # Check if pickle path exists, if it does exist, warn user that 
    # it will be overwritten
    if os.path.exists(output_pickle_path):
        print 'Pickle file %s already exists. \nPress enter to overwrite' \
            %output_pickle_path
        raw_input()
    
    # index uniprot database
    uniprot_dic, count = index_uniprot_database(uniprot_file, species='homosapiens')
    print '%s features indexed from %s' %(count, uniprot_file)
    
    # store dic as pickle
    pkl = open(output_pickle_path, 'wb')        
    pickle.dump(uniprot_dic, pkl)
    print 'Dictionary saved to: %s' %output_pickle_path
    
    
if __name__ == '__main__':
    main()