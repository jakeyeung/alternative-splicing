'''
Created on 2013-11-06

@author: jyeung
After running meme_batch_mode.sh, a large list of textfiles are created from MEME and TomTom.
Read the TomTom outputs matching RBP IDs and (after filtering) return RBP names.
'''

import sys
import csv
from optparse import OptionParser

def update_dic_with_rbp_info(mydic, rbp_row):
    '''
    Given a row in RBP_Information_all_motifs.txt (RBP DB file),
    get DBID as Key, and info as Values.
    Extract the following info:
    Example of DBID: RNCMPT00038.
    Example of DBID but not the one we want: ENSG00000076770
    Information to get (column names).
    Ensembl Gene ID, RBP_Name, RBP_Species, RBP_Status, Family_Name, RBDs. 
    
    Note: RBP DB file contains two columns with name DBID, therefore,
    we will use hard coded index column number to extract DBID.
    
    Justification for hard coding indexes rather than column names:
        There are duplicate column names, namingly DBID shows up twice.
    '''
    # Define my column name indices.
    # Columns containing my values.
    ensembl_index = 5
    rbp_name_index = 6
    rbp_species_index = 7
    rbp_status_index = 8
    family_name_index = 9
    rbd_index = 10
    # Columns containing my key (DBID)
    dbid_index = 11
    
    # Add DBID as key.
    if rbp_row[dbid_index] not in mydic:
        mydic[rbp_row[dbid_index]] = []    # create empty list.
    
    # Append values to list.
    for i in [ensembl_index, rbp_name_index, rbp_species_index, 
              rbp_status_index, family_name_index, rbd_index]:
        mydic[rbp_row[dbid_index]].append(rbp_row[i])
    return mydic

def index_rbp_file(rbp_db_path):
    '''
    Read RBP DB file, index based on DBID.
    Key in dic are DBID, values are relevant information
    (see update_dic_with_rbp_info for more info on values)
    '''
    # Init dictionary
    rbp_index_dic = {}
    with open(rbp_db_path, 'rb') as rbp_file:
        rbp_reader = csv.reader(rbp_file, delimiter='\t')
        rbp_reader.next()    # Skip header.
        for row in rbp_reader:
            rbp_index_dic = update_dic_with_rbp_info(rbp_index_dic, row)
    return rbp_index_dic

def main():
    usage = 'usage: %prog tomtom_file.txt rbp_database.txt output_file.txt'
    parser = OptionParser(usage=usage)
    (_, args) = parser.parse_args()
    if len(args) < 3:
        print('Three positional arguments required to '\
              'be specified in command line:\n'
              '1) TomTom output file, containing RBP IDs.\n'
              '2) RBP Database file, containing RBP information.\n'
              '3)Output file.')
        sys.exit()
    tomtom_path = args[0]
    rbp_db_path = args[1]
    output_path = args[2]
    
    # Index RBP DB file.
    rbp_index_dic = index_rbp_file(rbp_db_path)
    
    # 
    

if __name__ == '__main__':
    main()