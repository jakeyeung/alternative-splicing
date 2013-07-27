'''
Created on 2013-07-26

@author: jyeung

Reads the summary file outputted from bayes factor analysis and outputs
five regions of relevance:
    Upstream exon, skipped exon, downstream exon.
    Upstream intron and downstream intron.
Creating a bedfile of coordinates that can be loaded to UCSC browser and 
their sequences obtained through table browser.

This code has only been tested on skipped exon events...

The purpose of this code is to create coordinates to be loaded to UCSC browser,
their sequences obtained via UCSC browser to be sent for motif search.
'''


import sys
import csv
from read_write_util import read_write


def write_coords_to_multi_files(output_bed_file, bed_header, 
                                chrom, exon_starts, exon_ends, event,
                                score):
    '''
    List of exon-starts and ends, write to multiple files, one for each
    exon and junction. 
    '''
    # Check lengths
    if len(exon_starts) != len(exon_ends):
        print('Warning, length of exons starts and ends are not equal!'\
         '(%s vs. %s)' %(len(exon_starts), len(exon_ends)))
    
    # Initialize list of indices
    exon_index = range(0, len(exon_starts))
    junction_index = range(0, len(exon_starts) - 1)
    # Create keynames for exon and junction dictionaries.
    exon_keys = [''.join(['exon_', i]) for i in exon_index]
    junction_keys = [''.join(['juc_', i]) for i in junction_index]
    
    # Create dic to store output write objects to write exon coords
    output_file_dic_exons = {}
    for i, k in zip(exon_index, exon_keys):
        # Append output_bed_file suffix
        output_bed_file_suffixed = \
            ''.join([output_bed_file, '_exon', str(i)])
        output_file_dic_exons[k] = open(output_bed_file_suffixed, 'wb')
    
    # Create dic to store write objects to write junction coords
    output_file_dic_jucs = {}
    for i, k in zip(junction_index, junction_keys):
        # Append output_bed_file suffix
        output_bed_file_suffixed = \
            ''.join([output_bed_file, '_junction', str(i)])
        output_file_dic_jucs[k] = open(output_bed_file_suffixed, 'wb')
                   
    # Write exons to exon files.
    for start, stop, k in zip(exon_starts, exon_ends, exon_keys):
        # Create tab separated string for writing to file.
        write_str = '\t'.join([chrom, start, stop, event, score])
        output_file_dic_exons[k].write(write_str)
        
    # Write junctions to junc files.
    

def extract_coordinates_from_miso_bf(miso_file, output_bed_file, bed_header):
    '''
    Input:
        miso_file - file from bayes analysis. Expects column names such as
            event_name, chrom, and strand.
        bed_header - first line of bed file
    Output:
        output_bed_file - bed file extracting coordinates of upstream, skipped
        and downstream exon as well as the two upstream/downstream introns.
        
        Each region (upstream skipped, downstream, up/down introns) are outputted
        into a separate bed file, with the output_bed_file suffixed to indicate
        the type of region.
    '''
    # Define column name constants
    event_name_str = 'event_name'
    chrom_str = 'chrom'
    strand_str = 'strand'
    shading = '999'    # For UCSC
    
    # Create read_write_obj for reading and writing simultaneously and also works
    # in python 2.6
    read_write_obj = read_write(miso_file, output_bed_file, header=True)
    with open(miso_file, 'rb') as readfile:
        reader = csv.reader(readfile, delimiter='\t')
        # Get header
        colnames = reader.next()
        for row in read_write_obj.reader:
            event_name = row[colnames.index(event_name_str)]
            chrom = row[colnames.index(chrom_str)]
            strand = row[colnames.index(strand_str)]
            '''
            Event name needs to be parsed to get coordinates.
            
            Example event_name is: 
            "chr8:144679518:144679845:-@chr8:144671161:144672251:-
            @chr8:144668899:144669022:-"
            indicating upstream_exon:-@skipped_exon:-@downstream
            in the (-) strand.
            
            Strategy to extract coordinates:
            split by the colon (':') to get a list, 'l'.
            then determine number of exons by:
            # of exons = (l-1) / 3 because each exon
            has triplet of info (chrom, start, end) and 
            subtract one because last element is either (-) or (+).
            
            Count by groups of 3 along the list, 'l', extracting
            the second (start) and third (stop) element of each group.
            Ignore first element because it is chrom. 
            '''
            event_name_list = event_name.split(':')
            # Loop through event_name_list, get every 2nd and 3rd element
            # as a start/end pair.
            exon_starts = []
            exon_ends = []
            for i, coordinate, in event_name_list:
                if i % 3 == 0:
                    pass    # no remainder means you are at first element/chrom
                elif i % 3 == 1:
                    exon_starts.append(coordinate)
                elif i % 3 == 2:
                    exon_ends.append(coordinate)
                else:
                    sys.exit('Parsing error in %, %' %(i, coordinate))
            
            

def main():
    pass


if __name__ == '__main__':
    main()