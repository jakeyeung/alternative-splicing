'''
Created on 2013-09-06

@author: jyeung

Given a text file containing a column of events, append beside 
the column th gene names found from MISO annotations v2 hg19.
'''


import csv
import re
from optparse import OptionParser, Option


class write_obj(object):
    '''
    For reading and writing gene information.
    '''
    def __init__(self, output_path):
        '''
        Constructor
        '''
        self.writefile = open(output_path, 'wb')
        self.writeobj = csv.writer(self.writefile, delimiter='\t')
        self.path = output_path
    
    def close(self):
        self.writefile.close()

def create_csv_obj(opened_file, skip_header=True):
    '''
    Take opened_file like file = open(path, 'rb') and
    return its csv obj.
    '''
    reader = csv.reader(opened_file, delimiter='\t')
    if skip_header == True:
        reader.next()
        myheader = None
    else:
        myheader = reader.next()
    return reader, myheader

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
    
def get_event_gsymbol_from_row_str(row, ID_str, gsymbol_str, annot_dic):
    '''
    Grabs event and genesymbol then writes to annot_dic.
    '''
    event = get_info_from_annot_row(row, ID_str)
    genesymbol = get_info_from_annot_row(row, gsymbol_str)
    annot_dic = update_dic_with_info(annot_dic, event, genesymbol)
    return annot_dic

def get_event_gsymbol_from_row_index(row, ID_index, gsymbol_index, annot_dic):
    '''
    Grabs event and genesymbol then writes to annot_dic using index.
    '''
    '''
    first index should return something like
    ID=chr2:9624561:9624679:+@chr2:9627585:9627676:+@chr2:9628276:9628591:+
    '''
    event = get_info_from_annot_row(row, list_index=ID_index)
    '''
    list_index = -1 should return something like:
    gsymbol=IAH1
    '''
    genesymbol = get_info_from_annot_row(row, list_index=gsymbol_index)
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
    gsymbol_index = -1
    ID_index = 1
    gene_str = 'gene'
    annot_dic = {}
    rowcount = 0
    genecount = 0
    
    for annot_row in csv_readobj:
        if contains_mystring(annot_row, genestring_col_index, gene_str):
            annot_dic = get_event_gsymbol_from_row_index(annot_row, ID_index, gsymbol_index, annot_dic)
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
        myreader, _ = create_csv_obj(myfile, skip_header=True)
        annot_dic = get_dic_from_miso_reader(myreader)
    return annot_dic

def get_header(csv_read_obj):
    '''
    From csv_read_obj, read first row and store it as a header.
    '''
    return(csv_read_obj.next())

def match_event_to_gsymbol(event, annot_dic):
    '''
    Match an event to the gene symbol of annot_dic (key).
    '''
    try:
        gsymbol = annot_dic[event]
    except KeyError:
        gsymbol = None
    return gsymbol

def get_event_from_header(row, header):
    '''
    From row and given header, try to get event.
    It may either be 'event' or 'event_name'.
    No other possibilities!
    '''
    # define constants
    event_str = 'event'    # A column name in input_header
    event_str2 = 'event_name'    # Second try for event_name
    
    try:
        event = row[header.index(event_str)]
    except ValueError:
        event = row[header.index(event_str2)]
    return event

def iterate_rows_and_write_gsymbol(reader, header, annot_dic, writer_obj):
    '''
    Iterate rows in input_reader (a csv_read_obj) and try to
    find its corresponding gene name from dictionary.
    
    Searches for event_str in input_header, then uses that 
    as index for its iterations to extract an event from row.
    
    Writes new row to writefile.
    '''
    # Define constants
    original_header = header    # To preserve original index positions.
    header.insert(1, 'gsymbol')
    
    # Write inserted header to file.
    writer_obj.writerow(header)
    
    write_count = 0
    # Write header
    for row in reader:
        event = get_event_from_header(row, original_header)
        gsymbol = match_event_to_gsymbol(event, annot_dic)
        row.insert(1, gsymbol)
        writer_obj.writerow(row)
        write_count += 1
    return write_count

def main(input_filepath, annot_filepath, output_filepath):
    '''
    Strategy: read miso annotations, get dictionary of
    {event: genesymbol}, then run through input_filepath
    and append gene symbol for every event, output to new file
    called output_filepath.
    '''
    annot_dic = index_miso_annots(annot_filepath)
    
    # Initialize my writefile
    writer_obj = write_obj(output_filepath)
    
    with open(input_filepath, 'rb') as inputfile:
        input_reader, input_header = create_csv_obj(inputfile, 
                                                    skip_header=False)
        write_count = iterate_rows_and_write_gsymbol(input_reader, 
                                                     input_header, 
                                                     annot_dic, 
                                                     writer_obj.writeobj)
    print('%s rows written to file: %s' %(write_count, output_filepath))
    writer_obj.close()
    
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-i', '--input_filepath', dest='input_filepath',
                      help='Input file, likely the summary of t-test results'\
                      ' from a previous script.')
    parser.add_option('-a', '--annotation_filepath', dest='annot_filepath',
                      help='File of miso annotations (.gff3 file)')
    parser.add_option('-o', '--output_filepath', dest='output_filepath',
                      help='Output file name')
    (options, _) = parser.parse_args()
    
    main(options.input_filepath, options.annot_filepath, options.output_filepath)