'''
Created on 2013-11-14

@author: jyeung

Take a list of RBPs found in the RBPID database, then 
find the motifs of each RBP. Grab each motif and format them
so that they will be usable as input for FIMO.
'''

import sys
import os
import csv
from optparse import OptionParser
from utilities import writing_utils, reading_utils

def update_dic_with_row_info(row, header, rbp_db_dic):
    '''
    From row of RBP information all motifs file,
    get DBID, MOTIFID, RBPSTATUS (in order)
    
    rbpstatus can only be 'N' (no motif info), 
    'I' (indirect motif info) or 'D' (direct motif info)...
    
    We will only index rbps that have 'I' or 'D'...
    
    Final form looks like:
    {ensemblID: {RBP_Name:[], Motif_ID:[], RBP_Status:[]}}
    '''
    # Def constants
    rbpname_colname = 'RBP_Name'
    dbid_colname = 'DBID'
    motifid_colname = 'Motif_ID'
    rbpstatus_colname = 'RBP_Status'
    
    # Grab info from row using header and colname to grab index
    rbpname = row[header.index(rbpname_colname)]
    dbid = row[header.index(dbid_colname)]
    motifid = row[header.index(motifid_colname)]
    rbpstatus = row[header.index(rbpstatus_colname)]
    
    if rbpstatus == 'N':
        # No motif information, just return rbp_db_dic (do nothing).
        return rbp_db_dic
    
    elif rbpstatus == 'I' or rbpstatus == 'D':
        # Update RBP DB DIC using DBID as key, motifid and rbpstatus as subdic.
        if dbid not in rbp_db_dic:
            # Initialize key if key does not yet exist.
            rbp_db_dic[dbid] = {}
            # Add subkeys with empty lists as values.
            for subkey in rbpname_colname, motifid_colname, rbpstatus_colname:
                rbp_db_dic[dbid][subkey] = []
        # Add motifid and rbpstatus info to dic...
        for subkey, subvalue in zip([rbpname_colname, motifid_colname, rbpstatus_colname], 
                                    [rbpname, motifid, rbpstatus]):
            rbp_db_dic[dbid][subkey].append(subvalue)
        return rbp_db_dic
    
    else:
        print 'RBP_Status expected to be either "N", "I" or "D", %s found.' \
            %rbpstatus
        sys.exit()
    
def index_rbpdb_motifids(rbp_db_path):
    '''
    Read RBP DB file, grab the RBPID gene name (e.g. SRSF1)
    then 
    '''
    rbp_db_dic = {}
    with open(rbp_db_path, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        # Get first line as header.
        header = myreader.next()
        for row in myreader:
            rbp_db_dic = update_dic_with_row_info(row, header, rbp_db_dic)
    print '%s RBPs indexed.' %len(rbp_db_dic.keys())
    return rbp_db_dic

def get_motifids_from_rbp(ensemblid, rbp_db_dic, 
                          motif_subkey='Motif_ID',
                          rbpname_subkey='RBP_Name',
                          rbpstatus_subkey='RBP_Status'):
    '''
    RBP_DB_DIC is of form:
    {ensemblID: {RBP_Name:[], Motif_ID:[], RBP_Status:[]}}
    Match ensemblID, retrieve the motif id
    
    motif subkey is default "Motif_ID" for accessing the subkey to get motifids
    (same with rbpname and rbpstatus). 
    
    Return motif_list and motif_name, which will be defined as:
    RBPname,Motif_ID,RBP_Status
    '''
    try:
        motif_list = rbp_db_dic[ensemblid][motif_subkey]
    except KeyError:
        return [], []
    '''
    Define name of motif as RBPname,Motif_ID,RBP_Status
    First grab RBPname, motifIDs and RBPstatuses as lists
    Second, merge the lists into one list such that
    it is [RBPname1,Motif_ID1,Status1, RBPname2,ID2,Status2...]
    '''
    # 1.
    rbpname_list = rbp_db_dic[ensemblid][rbpname_subkey]
    rbpstatus_list = rbp_db_dic[ensemblid][rbpstatus_subkey]
    # 2.
    rbp_motif_name_list = []
    for rbpname, motifid, status in zip(rbpname_list, motif_list, rbpstatus_list):
        rbp_motif_name_list.append(','.join([rbpname, motifid, status]))
    return motif_list, rbp_motif_name_list

def main():
    usage = 'usage: %prog rbp_list.txtfile rbp_db_directory '\
        'meme_db_output.outputfile\n'\
        'Three args must be specified in commandline: \n'\
        '1) List of RBPs (ENSEMBL gene id).\n'\
        '2) RBP directory containing RBP_Information_all_motifs.txt'\
        ' and pwms_all_motifs directory.\n'\
        '3) Output file.\n'\
        '-h to display this help message.\n'
    parser = OptionParser(usage=usage)
    parser.add_option('-f', '--dbfile', dest='rbp_filename',
                      help='name of RBP DB file, '\
                      'default is RBP_Information_all_motifs.txt',
                      default='RBP_Information_all_motifs.txt')
    parser.add_option('-d', '--pwmdir', dest='pwm_dir',
                      help='Name of directory containing motif files.'\
                      ' Default is "pwms_all_motifs".',
                      default='pwms_all_motifs')
    (options, args) = parser.parse_args()
    if len(args) < 3:
        print 'Incorrect number of parameters specified.'
        print usage
        sys.exit()
    input_rbps_path = args[0]
    rbp_db_dir = args[1]
    output_path = args[2]
    rbp_filename = options.rbp_filename
    pwm_dir = options.pwm_dir
    
    # Check that writefile does not already exist.
    # if it does, then do not overwrite and just exit
    if os.path.isfile(output_path):
        print '%s already exists. Aborting custom MEME db creation.' \
            %output_path
        sys.exit() 
    
    # Get rbp_db_path
    rbp_db_path = os.path.join(rbp_db_dir, rbp_filename)
    
    # Read list of RBPs, store as a list.
    rbps_list = reading_utils.extract_column_from_textfile(input_rbps_path, 
                                                           col_to_extract=0)
    
    # Index RBP DB: {ensemblID: {RBP_Name:[], Motif_ID:[], RBP_Status:[]}}
    rbp_db_dic = index_rbpdb_motifids(rbp_db_path)
    
    # init write file
    with open(output_path, 'wb') as writefile:
        # Write headers
        writing_utils.write_meme_headers(writefile)
        # Match rbp ensemblID to indexed dic key, write motifs to file.
        motifcount = 0
        for rbp in rbps_list:
            motifid_list, motifname_list = get_motifids_from_rbp(rbp, 
                                                                 rbp_db_dic)
            for motifid, motifname in zip(motifid_list, motifname_list):
                # Add .txt suffix to motifid
                motifid_fname = ''.join([motifid, '.txt'])
                motifid_path = os.path.join(rbp_db_dir, pwm_dir, motifid_fname)
                # Read motifid path, get the motif.
                extracted_motif = \
                    reading_utils.read_motifs_from_file(motifid_path, 
                                                        skipheader=True,
                                                        rownames=True)
                # Only write something to file if extracted_motif is not empty.
                if len(extracted_motif) != 0:
                # Write extracted motif to file...
                    writing_utils.write_motif_to_file(writefile, 
                                                      extracted_motif, 
                                                      motifname, 
                                                      alength=4, 
                                                      nsites=20)
                    motifcount += 1
    print '%s motifs written to file: %s' %(motifcount, output_path)
    

if __name__ == '__main__':
    main()