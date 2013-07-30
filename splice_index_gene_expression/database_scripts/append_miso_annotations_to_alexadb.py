'''
Created on 2013-07-26

@author: jyeung
Using Jamie's existing alexaDB (WTASP file), append miso annotations.
We will focus on:
    SE
    MXE
    Alt3'SS
    Alt5'SS
'''


import sys
import csv

def write_jucs_to_dic(juc_dic, eventid, chromo, juc_start, juc_end, 
                      strand, juc_type):
    '''
    From junction start ends and other info, record it to dictionary
    But make sure no duplicates in entries. 
    '''
    # Write to dictionary, making sure no duplicate entries.
    juc_key = ':'.join([chromo, juc_start, juc_end])
    if juc_key not in juc_dic:
        juc_dic[juc_key] = {'eventid': [eventid], 
                                  'start': [juc_start], 
                                  'end': [juc_end], 
                                  'chromo': [chromo], 
                                  'strand': [strand],
                                  'type': [juc_type]}
    else:
        '''
        If exists already, append to list.
        '''
        for subkey, subval in zip\
            (['eventid', 'start', 'end', 'chromo', 'strand', 'type'],
             [eventid, juc_start, juc_end, chromo, strand, juc_type]):
            juc_dic[juc_key][subkey].append(subval)
    return juc_dic

def get_coords_from_eventid(juc_coord_dic, eventid, chromo, strand, 
                            eventtype='SE'):
    '''
    From MISO annotations, get exon coordinates from eventid.
    # Parse event id to get upstream, skipped exon 
    # and downstream exons.
    
    Output as junction coordinates in a dictionary.
    
    Dictionary has the following keynames:
    eventid, start, end, chromo, strand, type
    '''
    # Initialize constants
    
    # Take first element in semi-colon separated string.
    # Then remove first three elements (contains 'ID=')
    eventid_parsed_id = eventid.split(';')[0][3:]    
    # Split by colon and remove last element in list
    # because it is strand, which we already know.
    eventid_parsed = eventid_parsed_id.split(':')[:-1]
    # Loop through eventid, get the strings that do not contain chromo.
    coords = []
    for i in eventid_parsed:
        if chromo not in i:
            '''
            Order of coords begins from start of transcript to end.
            But coord pairs (exon start, exon end) is such that 
            exon end > exon start.
            
            Therefore, for negative strand, it looks like:
            1000, 1500, 500, 750, 100, 250
            
            Whereas for positive strand it looks like:
            100, 250, 500, 750, 1000, 1500
            
            So be careful!
            '''
            coords.append(i)
            
    #raw_input(coords)
    
    # Should be even number.
    if len(coords) % 2 != 0:
        print('Warning, non-even number of coords, cannot pair up.')
        sys.exit()
    # BEGIN: Get junction starts and ends
    # Calculate number of exons
    n_exons = len(coords) / 2
    n_jucs = n_exons - 1
    for i in range(0, n_jucs):
        # Get juc start and ends, depending on strand
        juc_count = 0    # for SE events.
        if strand == '+':
            juc_start = coords[(i * 2) + 1]
            for j in range(i+1, n_exons):
                juc_end = coords[(j * 2)]
                if eventtype == 'SE':
                    if juc_count == 0:
                        juc_type = 'inclusion'
                    elif juc_count == 1:
                        juc_type = 'exclusion'
                    else:
                        print('For SE events, juc_count, %s, '\
                              'not expected to be >1.' %juc_count)
                        sys.exit() 
                juc_coord_dic = \
                    write_jucs_to_dic(juc_coord_dic,        
                                      eventid_parsed_id, chromo, 
                                      juc_start, juc_end, 
                                      strand, juc_type)
                '''
                print(strand)
                print(i, j)
                print((juc_start, juc_end))
                raw_input(coords)
                '''
                juc_count += 1
        elif strand == '-':
            #print range(0, n_jucs)
            #raw_input(range(i+1, n_jucs)) 
            juc_end = coords[(i * 2)]
            for j in range(i+1, n_exons):
                juc_start = coords[(j * 2) + 1]
                if eventtype == 'SE':
                    if juc_count == 0:
                        juc_type = 'inclusion'
                    elif juc_count == 1:
                        juc_type = 'exclusion'
                    else:
                        print('For SE events, juc_count, %s, '\
                              'not expected to be >1.' %juc_count)
                        sys.exit()
                juc_coord_dic = \
                    write_jucs_to_dic(juc_coord_dic, 
                                      eventid_parsed_id, chromo, 
                                      juc_start, juc_end, 
                                      strand, juc_type)
                print ':'.join([chromo, juc_start, juc_end])
                print juc_coord_dic[':'.join([chromo, juc_start, juc_end])]
                print(strand)
                print(i, j)
                print((juc_start, juc_end))
                raw_input(coords)
                juc_count += 1
        
                
    return juc_coord_dic

def index_annotations(annot_file, eventtype='SE'):
    '''
    Open miso annotation file and create indexed dictionaries.
    
    Dictionary format: 
        key: colon joined chromo, juc_start, juc_end
        value: dictionary with keys: eventid, start, end, chromo, strand, type 
    '''
    # Initialize constants
    chr_index = 0
    eventtype_index = 1
    gene_or_exon_index = 2
    strand_index = 6
    id_index = 8
    
    # Initialize dics
    juc_dic = {}
    
    if eventtype == 'SE':
        '''
        For skipped events:
        Exon 1 means upstream constitutive exon.
        Exon 2 means skipped exon.
        Exon 3 means downstream constitutive exon.
        
        In Exon 1, we will record the location of exon 1's stop coordinate.
        In Exon 2, we will record both start and stop coordinates.
        in Exon 3, we will record only exon's start coordinate.
        
        Rationale:
            We will search for junctions in alexaDB. These junctions 
            will only connect in three ways:
                Exon 1 to Exon 2
                Exon 2 to Exon 3
                Exon 1 to Exon 3.
                
        We also need to be wary of negative and positive strands, which
        may reverse the order of start and stop for exons 1, 2, 3.
        '''
    rowcount = 0
    with open(annot_file, 'rb') as readfile:
        reader = csv.reader(readfile, delimiter='\t')
        # Ignore header
        reader.next()
        '''
        MISO annotation format:
        chr1, SE, gene/mRNA/exon, start, stop, ., strand, ., ID
        Where ID is:
        ID=chr1:4775654:4775821:-@chr1:4774032:4774186:-@chr1:4772649:4772814:-;
        Name=chr1:4775654:4775821:-@chr1:4774032:4774186:-@chr1:4772649:4772814:-
        '''
        for row in reader:
            # Only look at row if third column == 'gene'
            if row[gene_or_exon_index] != 'gene':
                pass
            else:
                # Get values of the gene that shows differential expression
                if row[eventtype_index] == eventtype:
                    pass
                else:
                    print('Warning, expected event type in annotation '\
                          'file to be: %s' %eventtype)
                    print('Continuing anyway...')
                chromo = row[chr_index]
                eventid = row[id_index]
                strand = row[strand_index]
                # Get start and end coordinates for junctions.
                juc_dic = \
                    get_coords_from_eventid(juc_dic, eventid, chromo, strand)
                rowcount += 1
    return juc_dic, rowcount
                
                


def main():
    if len(sys.argv) < 5:
        print('Some parameters need to be specified in command line.')
        sys.exit()
    alexa_bed_file = sys.argv[1]
    annot_file = sys.argv[2]    # From MISO
    output_bed_file = sys.argv[3]
    eventtype = sys.argv[4]
    
    juc_dic, rowcount = \
        index_annotations(annot_file, eventtype=eventtype)
    
    # print juc_dic.keys()
    print len(juc_dic.keys())
    
    

if __name__ == '__main__':
    main()
    