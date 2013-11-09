'''
Created on 2013-11-06

@author: jyeung
After running meme_batch_mode.sh, a large list of textfiles are created from MEME and TomTom.
Read the TomTom outputs matching RBP IDs and (after filtering) return RBP names.
'''

import sys
import csv
from optparse import OptionParser
from utilities import writing_utils

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
    dbid_index = 13
    
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

def passes_orientation(row, header, filter_strand):
    '''
    If filter_strand is True:
        Given row, checks if row orientation is + (boolean)
    Otherwise, returns True no matter what.
    '''
    # Def colname
    orientation_colname = 'Orientation'
    if filter_strand:
        return row[header.index(orientation_colname)] == '+'
    else:
        return True

def passes_qval_cutoff(row, header, filter_qval):
    '''
    Given row, checks if qvalue is less than or equal to
    filter_qval.
    '''
    # def colname
    qval_colname = 'q-value'
    qval = float(row[header.index(qval_colname)])
    return qval <= filter_qval

def get_rbp_info(llists, row, header, rbp_index_dic):
    '''
    Get relevant RBP information appending it as a list to
    llist (making a list of lists).
    
    First: get info from row.
    Second: get info from rbp_index_dic
    '''
    # Def colnames
    # TomTom colnames 
    targetid_colname = 'Target ID'
    qval_colname = 'q-value'
    overlap_colname = 'Overlap'
    query_colname = 'Query consensus'
    target_colname = 'Target consensus'
    
    # Begin getting info from RBP dictionary.
    rbp_info = []
    rbp_id_key = row[header.index(targetid_colname)]
    try:
        rbp_info += rbp_index_dic[rbp_id_key]
        # Begin getting info from row.
        for col in [targetid_colname, qval_colname, overlap_colname, 
                    query_colname, target_colname]:
            rbp_info.append(row[header.index(col)])
        llists.append(rbp_info)
    except KeyError:
        # print '%s not in human RBP database.' %rbp_id_key
        pass
    return llists

def extract_rbps_from_tomtom(tomtom_path, rbp_index_dic, 
                             filter_qval=1, 
                             filter_strand=True):
    '''
    Reads tomtom path, matches TargetID from TomTom to rbp dictionary.
    Extract:
        1) TargetID
        2) Q value
    From RBP ID extract values corresponding to TargetID.
    Options:
        filter_qval: only looks at hits with qval less
             than or equal to filter_qval
        filter_strand: if true: only look at positive strand orientations.        
    '''
    rbp_matches = []    # list of lists.
    with open(tomtom_path, 'rb') as tomtom_file:
        tomtom_reader = csv.reader(tomtom_file, delimiter='\t')
        tomtom_header = tomtom_reader.next()
        for row in tomtom_reader:
            if passes_orientation(row, tomtom_header, filter_strand) and \
                passes_qval_cutoff(row, tomtom_header, filter_qval):
                # Get info from TomTom file and indexed dictionary.
                rbp_matches = get_rbp_info(rbp_matches, row, tomtom_header, 
                                           rbp_index_dic)
    return rbp_matches

def main():
    usage = 'usage: %prog tomtom_file.txt rbp_database.txt output_file.txt'
    parser = OptionParser(usage=usage)
    parser.add_option('-q', '--q_cutoff', dest='qval_cutoff', default=0.15,
                      help='Specify Q-value cutoff, default 0.15.')
    parser.add_option('-s', '--filter_strand', dest='filter_strand', 
                      default='True',
                      help='Specify whether to filter strand (True->"+" only)')
    (options, args) = parser.parse_args()
    if len(args) < 3:
        print('Three positional arguments required to '\
              'be specified in command line:\n'
              '1) TomTom output file, containing RBP IDs.\n'
              '2) RBP Database file, containing RBP information.\n'
              '3) Output file.')
        sys.exit()
    tomtom_path = args[0]
    rbp_db_path = args[1]
    output_path = args[2]
    qval_cutoff = float(options.qval_cutoff)
    if options.filter_strand in ['True', 'true', 'TRUE', 'T', 't']:
        filter_strand = True
    elif options.filter_strand in ['False', 'false', 'FALSE', 'F', 'f']:
        filter_strand = False
    else:
        print 'Filter strand must be True or False. %s found.' \
        %options.filter_strand
        sys.exit()
     
    # Index RBP DB file.
    rbp_index_dic = index_rbp_file(rbp_db_path)
    
    # Read TomTom file, extract relevant RBPs.
    matched_rbps = extract_rbps_from_tomtom(tomtom_path, 
                                            rbp_index_dic,
                                            filter_qval=qval_cutoff,
                                            filter_strand=filter_strand)
    # Write TomTom file to output.
    writing_utils.write_list_of_lists_to_file(matched_rbps, output_path)

if __name__ == '__main__':
    main()