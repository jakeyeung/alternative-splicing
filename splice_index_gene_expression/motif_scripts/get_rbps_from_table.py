'''
Created on 2013-10-20

@author: jyeung

Downloaded RBP Motifs from cisrna-bp database
http://cisbp-rna.ccbr.utoronto.ca/

I got a table of RBP information for homo sapiens: 619 rows.

Parse this file and grab 
'''

import csv
import sys
import os
from optparse import OptionParser


def loop_rows_append_rbp_dic(rbp_dic, myreader):
    '''
    Input: 
    rbp_dic: should contain empty lists. {D:[], I:[], N:[]}
    myreader: csv object of rbp motif information from cisbp
    
    Add ensembl ID followed by RBI name to the dictionary.
    '''
    # Initialize constants
    rbp_status_colname = 'RBP_Status'
    rbp_name_colname = 'RBP_Name'
    ensembl_id_colname = 'DBID'
    
    header = myreader.next()    # First row is colnames.
    for row in myreader:
        rbp_name = row[header.index(rbp_name_colname)]
        rbp_id = row[header.index(ensembl_id_colname)]
        rbp_status = row[header.index(rbp_status_colname)]
        try:
            rbp_dic[rbp_status].append((rbp_id, rbp_name))
        except KeyError:
            print('Unknown RBP status. '\
            'Expected "D", "I", or "N". %s found.' %rbp_status)
    # Get set to remove redundant genes.
    for k in rbp_dic.keys():
        rbp_dic[k] = list(set(rbp_dic[k]))
    return rbp_dic
    
def index_rbps_by_rbp_status(input_path):
    '''
    Loop through RBP motif data, output a dictionary
    with keys as rbp status. 3 rbp status possible:
    Where D: directly determined motif. 
          I: inferred from another RBP, based on RBP similarity.
          N: no motif is available
    Example output:
      {D:[rbpsD], I:[rbpsI], N:[rbpsN]}
    '''
    
    # Initialize output dic with empty lists.
    rbp_dic = {'D':[], 'I':[], 'N':[]}
    
    # Initialize read object then loop through rows.
    with open(input_path, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        rbp_dic = loop_rows_append_rbp_dic(rbp_dic, myreader)
        
    return rbp_dic

def add_suffix_to_path(path, suffix):
    '''
    Given a path, add a suffix to it.
    e.g.
    /Data/jyeung/myfile.txt -> /Data/jyeung/myfile_suffix.txt
    '''
    path_split = os.path.splitext(path)
    suffixed_path = ''.join([path_split[0], '_', suffix, path_split[1]])
    return suffixed_path

def write_list_to_file(mylist, mypath):
    '''
    Given a list, write to file.
    '''
    # Init my file, then loop list.
    writecount = 0    # number of rows written to file.
    with open(mypath, 'wb') as myfile:
        mywriter = csv.writer(myfile, delimiter='\t')
        for l in mylist:
            mywriter.writerow(l)
            writecount += 1
    return writecount
        

def write_rbps_to_file(rbp_dic, output_path):
    '''
    Will write three text files, one for Direct, Indirect, and None.
    Each text file will contain a list of RBPs.
    '''
    output_path_D = add_suffix_to_path(output_path, 'direct')
    output_path_I = add_suffix_to_path(output_path, 'indirect')
    output_path_N = add_suffix_to_path(output_path, 'not_in_db')
    
    for path, status in zip([output_path_D, output_path_I, output_path_N], 
                            ['D', 'I', 'N']):
        writecount = write_list_to_file(rbp_dic[status], path)
        print '%s rows written to: %s' %(writecount, path)
    return None

def main():
    usage = 'usage: %prog input_filepath output_path'
    parser = OptionParser(usage=usage)
    (_, args) = parser.parse_args()    # no options at the moment.
    
    if len(args) < 2:
        print('Input path and output path must be provided in command line. '\
              '\nType -h for help')
        sys.exit()
    input_path = args[0]
    output_path = args[1]
    
    rbp_dic = index_rbps_by_rbp_status(input_path)
    print('Direct RBPs: %s\nIndirect RBPs: %s\nN RBPs: %s'\
           %(len(rbp_dic['D']), len(rbp_dic['I']), len(rbp_dic['N'])))
    
    write_rbps_to_file(rbp_dic, output_path)
    
if __name__ == '__main__':
    main()