'''
Created on 2013-08-01

@author: jyeung

Retrieve gene names from GTF file given the exon coordinates. 
'''


import sys
import csv
import re

def search_between_quotes(text):
    '''
    Given a text, uses regex expressions to search whatever is between two 
    quotations, " "
    '''
    pattern = r'"([A-Za-z0-9_\./\\-]*)"'
    m = re.search(pattern, text)
    try:
        my_match = m.group()
    except AttributeError:
        print('Could not find anything between quotes.')
        print('Returning none...')
        my_match = None
    return my_match
    
def index_gtf_file(gtf_file):
    '''
    gtf file format (in GRCh37.72 ensembl):
        col1: chromosome (no 'chr' prefix)
        col4: start
        col5: end
        col7: strand
        col9: attribute info (semi-colon separated values: 
            example: gene_id; transcript_id; exon_number; gene_name;...)
    
    Create dictionaries of the form: 
        key: "chrX:start:end" <- searchable via miso id
        value: ensemblid, gene_name
    '''    
    # Define column index constants
    chromo_index = 0
    start_index = 3
    stop_index = 4
    strand_index = 6
    attrib_index = 8
    prefix = 'chr'    # For appending to ensembl chr name.
    
    # Define ensembl chromosome list we care about.
    chromo_list = [str(i) for i in range(1, 23) + ['X', 'Y']]
    
    exon_dic = {}
    with open(gtf_file, 'rb') as readfile:
        reader = csv.reader(readfile, delimiter='\t')
        for row in reader:
            chromo = row[chromo_index]
            # Only care about chromosomes 1, 2, 3..., 22, X, Y
            if chromo not in chromo_list:
                # Skip to next row
                continue
            chromo_prefixed = ''.join([prefix, chromo])
            start = row[start_index]
            stop = row[stop_index]
            strand = row[strand_index]
            # Create dic key like chr1:123:345:-
            dic_key = ':'.join([chromo_prefixed, start, stop, strand])
            # Check if dic_key already exists, if it does not
            # exist, initialize a list.
            if dic_key not in exon_dic:
                exon_dic[dic_key] = []
            # Get ensemblID and gene name by parsing semicolon split attrib
            attrib = row[attrib_index].split(';')
            # retrieve 1st and 4th element 
            # in split list (gene_id and gene_name, respectively).
            ensembl_id_str = attrib[0]
            gene_name_str = attrib[3]
            # Search between quotations " " to get id and name.
            ensembl_id = search_between_quotes(ensembl_id_str)
            gene_name = search_between_quotes(gene_name_str)
            # Append to list as a tuple (ensembl_id, gene_name)
            exon_dic[dic_key].append((ensembl_id, gene_name))
    return exon_dic

def read_miso_id_write_gene_name(miso_file, exon_dic, output_file):
    '''
    miso_ids can come from a source such as miso filtered bf output, which
    contains the column id "event_name". We will search this name to get id.
    
    exon_dic is output of index_gtf_file(), with form:
        key-> chr1:123:345
        value-> [(ensembl_id1, gene_name1), (ensembl_id2, gene_name2),...]
    
    Write all info from miso_file, but include last column to be gene name.
    '''
    # Initialize colname constants
    eventid_str = 'event_name'    # Should match colname in miso_file
    genename_str = 'gene_name'    # Colname for output file.
    exon_index = 1    # Which exon to take from eventid?
    
    # Open output file
    writefile = open(output_file, 'wb')
    writer = csv.writer(writefile, delimiter='\t')
    # Open miso file
    with open(miso_file, 'rb') as readfile:
        reader = csv.reader(readfile, delimiter='\t')
        header = reader.next()
        # Add gene name to end of header and write as first line.
        header.append(genename_str)
        writer.writerow(header)
        for row in reader:
            # Get event id, parsed by @
            eventid_parsed = row[header.index(eventid_str)].split('@')
            # In SE events, there will be 3 exons, take the middle one. 
            dic_key = eventid_parsed[exon_index]
            # Search dic key
            if dic_key in exon_dic:
                # tuple of form: (ensembl_id, gene_name)
                gene_name_list = [i[1] for i in list(set(exon_dic[dic_key]))]
                # Record as CSV
                row.append(','.join(gene_name_list))
                writer.writerow(row)
            else:
                pass
    writefile.close()

def main():
    if len(sys.argv) < 4:
        print('GTF file, miso_bf_filtered file and output file '\
              'must be specified in command line.')
        sys.exit()
    gtf_file = sys.argv[1]
    miso_file = sys.argv[2]
    output_file = sys.argv[3]
    
    # Index gtf file
    exon_dic = index_gtf_file(gtf_file)
    
    # Read miso file, get gene name, write to output.
    read_miso_id_write_gene_name(miso_file, exon_dic, output_file)
    
if __name__ == '__main__':
    main()