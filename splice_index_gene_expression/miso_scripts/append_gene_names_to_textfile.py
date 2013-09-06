'''
Created on 2013-09-06

@author: jyeung

Given a text file containing a column of events, append beside 
the column th gene names found from MISO annotations v2 hg19.
'''


import csv
import sys
import re


def create_csv_obj(opened_file, skip_header=True):
    '''
    Take opened_file like file = open(path, 'rb') and
    return its csv obj.
    '''
    reader = csv.reader(opened_file, delimiter='\t')
    if skip_header == True:
        reader.next()
    return reader

def create_regexp_from_searchstr(searchstr):
    '''
    Create reg expression such that we will search for
    everything after searchstr+.
    Example, if searchstr = 'ID', we want to search for
    everything after 'ID+'
    
    S+ means any non-white space character.
    '''
    regexp = ''.join(['(?<=', searchstr, '=)\S+'])
    # Example: regexp = (?<=ID=)\S+
    return regexp

def iterate_search_from_split_str_list(split_str_list, regexp):
    '''
    From split_str_list, try to find which one returns
    something that is NoneType.
    '''
    result = False
    for jstr in split_str_list:
        match = re.search(regexp, jstr)
        if match:
            result = match.group(0)
            break
    return result

def search_from_split_str(split_str_list, searchstr):
    '''
    Example: if searchstr = 'ID', we will get
    chr2:9624561:9624679:+@chr2:9627585:9627676:+@chr2:9628276:9628591:+
    from
    ID=chr2:9624561:9624679:+@chr2:9627585:9627676:+@chr2:9628276:9628591:+
    
    Note split_str_list is actually a LIST, basically iterate 
    '''
    regexp = create_regexp_from_searchstr(searchstr)
    result = iterate_search_from_split_str_list(split_str_list, regexp)
    if result == False:
        print('Warning: no match found for %s' %searchstr)
    return result

def search_from_split_str2(split_str_list, list_index):
    '''
    Instead of using regex, may be faster to just use list_index then
    split with '='.
    '''
    mystr = split_str_list[list_index].split('=')[1]
    return mystr

def get_info_from_annot_row(annot_row, list_index):
    '''
    From a row in miso annot file, get its event.
    Assumes a semicolon separated string containing event 
    in the 9th column (index 8).
    '''
    # Define constants
    annot_str_index = 8
    annot_str = annot_row[annot_str_index]
    # info = search_from_split_str(annot_str.split(';'), searchstr)
    info = search_from_split_str2(annot_str.split(';'), list_index)
    return info

def update_dic_with_info(mydic, mykey, myval):
    '''
    Create key:value pair but give warning if the key already exists
    in the dictionary.
    '''
    if mykey not in mydic:
        mydic[mykey] = myval
    else:
        print('Warning, %s already exists in dictionary, '\
              'overwriting existing values...' %mykey)
        mydic[mykey] = myval
    return mydic

def contains_mystring(annot_row, col_index, mystring):
    '''
    Checks ith column (e.g. 2) if string contains mystring (e.g. 'gene')
    Returns boolean True or False.
    '''
    return(annot_row[col_index] == mystring)
    
def get_event_gsymbol_from_row(row, ID_str, gsymbol_str, annot_dic):
    '''
    Grabs event and genesymbol then writes to annot_dic.
    '''
    '''
    first index should return something like
    ID=chr2:9624561:9624679:+@chr2:9627585:9627676:+@chr2:9628276:9628591:+
    '''
    event = get_info_from_annot_row(row, list_index=0)
    '''
    list_index = -1 should return something like:
    gsymbol=IAH1
    '''
    genesymbol = get_info_from_annot_row(row, list_index=-1)
    annot_dic = update_dic_with_info(annot_dic, event, genesymbol)
    return annot_dic

def get_dic_from_miso_reader(csv_readobj, only_gene=True):
    '''
    From csv_readobj, likely created from create_csv_obj,
    iterate its rows, collecting its eventname and genesymbol.
    
    only_gene = True: checks if third column has the string 'gene', 
    skips if it is not 'gene'. This is useful because
    gene symbol in miso annotation only is in 'gene', not 'mRNA'
    or 'exon'.
    '''
    # Define constants and empty dics
    genestring_col_index = 2
    gene_str = 'gene'
    ID_str = 'ID'
    gsymbol_str = 'gsymbol'
    annot_dic = {}
    rowcount = 0
    genecount = 0
    
    for annot_row in csv_readobj:
        if contains_mystring(annot_row, genestring_col_index, gene_str):
            annot_dic = get_event_gsymbol_from_row(annot_row, ID_str, 
                                                   gsymbol_str, annot_dic)
            genecount += 1
        else:    # not 'gene', then probably doesnt contain gsymbol.
            pass
        rowcount += 1
    
    print ('Iterated %s rows, of them %s were generows.' %(rowcount, genecount))
    return annot_dic
        
def index_miso_annots(annot_filepath):
    '''
    Opens file using 'with', creates csv_read_obj, 
    then runs functions to get annot_dic.
    '''
    with open(annot_filepath, 'rb') as myfile:
        myreader = create_csv_obj(myfile, skip_header=True)
        annot_dic = get_dic_from_miso_reader(myreader)
    return annot_dic

def main(input_filepath, annot_filepath, output_filepath):
    '''
    Strategy: read miso annotations, get dictionary of
    {event: genesymbol}, then run through input_filepath
    and append gene symbol for every event, output to new file
    called output_filepath.
    '''
    miso_dic = index_miso_annots(annot_filepath)
    
if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Textfile, miso annotations and output file must be specified '\
              'in command line.')
        sys.exit()
    input_filepath = sys.argv[1]
    annot_filepath = sys.argv[2]
    output_filepath = sys.argv[3]
    main(input_filepath, annot_filepath, output_filepath)