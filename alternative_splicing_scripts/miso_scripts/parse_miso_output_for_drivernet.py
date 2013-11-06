'''
Created on 2013-09-24

@author: jyeung

For drivernet, we require a table containing:
SampleX    Gene

For each differentially spliced event.

We can do this by parsing miso output.
'''

import sys
import csv

class csv_obj():
    '''
    Prepare read obj for reading miso output
    TODO: this is now in utilities, so you can delete this and use
    the function in utilities.
    '''
    def __init__(self, miso_output, filetype):
        self.filepath = miso_output
        self.type = filetype
        if self.type != 'read' and self.type != 'write':
            print('filetype must be either "read" or'\
                  ' "write". %s found...' %filetype)
            sys.exit()
        
    def __enter__(self):
        if self.type == 'read':
            self.file = open(self.filepath, 'rb')
            self.readobj = csv.reader(self.file, delimiter='\t')
        elif self.type == 'write':
            self.file = open(self.filepath, 'wb')
            self.writeobj = csv.writer(self.file, delimiter='\t')
        else:
            print('filetype must be either '\
                  '"read" or "write".. %s found.' %self.type)
            sys.exit()
    
    def __exit__(self, jtype, jvalue, tb):
        self.file.close()
        
def get_key_value_list_from_row(row, header, key_str, value_str):
    '''
    
    '''
    jkey = row[header.index(key_str)]
    jvalue = row[header.index(value_str)]
    jvalue_list = jvalue.split(',')
    return jkey, jvalue_list
        
def create_dic_from_readobj(key_str, value_str, header, readobj):
    '''
    Iterate rows from readobj, extract info from each row by 
    using key_str, value_str and header to get proper index.
    Returns a dic with mykey: myvalues_list.
    
    Expects myvalues to be in CSV form, so try to join to list.
    
    Example:
    {MAP4K4: [tophat_97_T_PE,tophat_2682_A_PE,tophat_2741_B_PE,
    tophat_2743_D_PE,tophat_3035_B53_PE,tophat_3071_B51_PE,
    tophat_7800_T76_PE,tophat_7820_PE]}
    
    If key already exists, then simply append to existing list.
    '''
    jdic = {}
    for row in readobj:
        jkey, jvalue_list = get_key_value_list_from_row(row, header, 
                                                        key_str, value_str)
        if jkey not in jdic:
            jdic[jkey] = jvalue_list
        else:    # Key already exists in dic
            jdic[jkey] += jvalue_list
            # Remove repeats from list...
            jdic[jkey] = list(set(jdic[jkey]))
    print('%s genes loaded onto dictionary.' %len(jdic.keys()))
    return jdic
        
def make_gene_sample_dic(readobj):
    '''
    creates read obj and returns a dictionary with
    {gene: samples_list} form. 
    
    readobj is a csv readobject, with headers as first row.
    '''
    # Def constants to match header names.
    gsymbol_str = 'gsymbol'
    sample_name_str = 'sample_name'
    # Create dictionary
    header = readobj.next()
    gene_sample_dic = create_dic_from_readobj(gsymbol_str, sample_name_str, 
                                              header, readobj)
    return gene_sample_dic

def write_header_to_file(header, writeobj):
    '''
    Write header to file
    '''
    writeobj.writerow(header)

def write_key_val_to_file(dic, writeobj):
    '''
    Write dic to file through writeobj.
    
    Since values are in a list form, we want each element in list
    to be in its own row.
    
    Example:
    val1: key1
    val2: key1
    val3: key1
    val1: key2
    etc.
    '''
    writecount = 0
    for key, value_list in dic.iteritems():
        for value in value_list:
            writeobj.writerow([value, key])
            writecount += 1
    return writecount

def main(miso_output, aberrant_output):
    miso_readobj = csv_obj(miso_output,filetype='read')
    with miso_readobj:
        gene_samp_dic = make_gene_sample_dic(miso_readobj.readobj)
        writeobj = csv_obj(aberrant_output, filetype='write')
        with writeobj:
            write_header_to_file(['sample', 'gene'], writeobj.writeobj)
            writecount = write_key_val_to_file(gene_samp_dic, 
                                               writeobj.writeobj)
    print('%s rows written to file: %s' %(writecount, writeobj.filepath))
            
if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('MISO output and output file must be specified in command line.')
        sys.exit()
    miso_output = sys.argv[1]
    aberrant_output = sys.argv[2]
    main(miso_output, aberrant_output)
    