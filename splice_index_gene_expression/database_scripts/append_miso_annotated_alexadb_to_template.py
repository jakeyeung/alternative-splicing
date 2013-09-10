'''
Created on 2013-09-10

@author: jyeung

Read output from append_miso_annotations_to_alexadb.py text file and store 
it as a dictionary.

Read template file and write new file containing the indexed information
'''


import sys
import csv

def add_values_to_subdic(index_dic, juc_id, myvalue, colname):
    '''
    Put values into subdic!
    '''
    if juc_id not in index_dic:    # intialize no key found.
        index_dic[juc_id] = {}
        index_dic[juc_id][colname] = myvalue
    else:
        index_dic[juc_id][colname] = myvalue
    return index_dic

def update_dic_with_row_info(row, header, index_dic, colname):
    '''
    grab row info using header to find row index.
    Update index_dic accordingly.
    
    if colname==eventid searches index of the row to obtain eventid
    adds it to index_dic.
    '''
    # Define constants
    junction_id_str = 'junction_id'
    juc_id = row[header.index(junction_id_str)]
    myvalue = row[header.index(colname)]
    index_dic = add_values_to_subdic(index_dic, juc_id, myvalue, colname)
    return index_dic

def store_rows_to_dic(readobj):
    '''
    Itereate rows in readobj, get information into dic.
    
    Assumes first row in readobj is a header.
    '''
    index_dic = {}    # initialize
    
    header = readobj.next()
    rowcount = 0
    for row in readobj:
        index_dic = update_dic_with_row_info(row, header, index_dic, 'eventid')
        index_dic = update_dic_with_row_info(row, header, index_dic, 'type')
        rowcount += 1
    print('%s junctions indexed to dictionary.' %rowcount)
    return index_dic

def index_miso_alexa(miso_alexa_bed):
    '''
    Read miso_alexa_bed file, store:
        -jucid (key)
        -eventid (subkey)
        -strand (subkey)
        -type (subkey)
    '''
    with open(miso_alexa_bed, 'rb') as myfile:
        myreader = csv.reader(myfile, delimiter='\t')
        index_dic = store_rows_to_dic(myreader)
    return index_dic

def initialize_output_file(output_file):
    '''
    Create write object
    '''
    outfile = open(output_file, 'wb')
    out_writer = csv.writer(outfile, delimiter='\t')
    return out_writer
    
def expand_header(header, new_columns):
    '''
    take header (list) append list to new columns. 
    '''
    return header + new_columns

def row_match_miso_juc(row, juc, index_dic):
    '''
    Check if row matches with any indexed junction
    '''
    if juc in index_dic:
        row_match = True
    else:
        row_match = False
    return row_match

def append_template_with_miso_info(index_dic, template_file, output_file):
    '''
    Strategy: read template file, iterate through rows and outputting each
    row to output file, but modifying the output to include index_dic
    if it matches.
    '''
    writecount = 0
    # Prepare read file
    with open(template_file, 'rb') as readfile:
        # Make read obj, get its header.
        reader = csv.reader(readfile, delimiter='\t')
        read_header = reader.next()
        # Prepare output file, write expanded header.
        out_writer = initialize_output_file(output_file)
        write_header = expand_header(read_header, ['eventid', 'type'])
        out_writer.writerow(write_header)
        # Iterate rows write to file.
        for row in reader:
            juc = row[read_header.index('id')]
            if row_match_miso_juc(row, juc, index_dic)==True:
                for subkey in ['eventid', 'type']:
                    subval = index_dic[juc][subkey]
                    row.append(subval)
                out_writer.writerow(row)
            else:    # No match, expand_row with None.
                out_writer.writerow(row + [None, None])
            writecount += 1
            if writecount % 10000==0:
                print('Through %s junctions...')
    print('%s rows written to file %s:' %(writecount, output_file))
    return writecount

def main(miso_alexa_bed, template_file, output_file):
    index_dic = index_miso_alexa(miso_alexa_bed)
    append_template_with_miso_info(index_dic, template_file, output_file)

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('miso-alexadb-bed file, template file and output file '
              'must be specified in command line.')
        sys.exit()
    miso_alexa_bed = sys.argv[1]    # from append_miso_annotations_to_alexadb.py
    template_file = sys.argv[2]
    output_file = sys.argv[3]
    main(miso_alexa_bed, template_file, output_file)