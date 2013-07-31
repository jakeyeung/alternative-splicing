'''
Created on 2013-07-26

@author: jyeung
Using Jamie's existing alexaDB (WTASP file), append miso annotations.
We will focus on:
    SE
    MXE
    Alt3'SS
    Alt5'SS
    
#TODO: Remove those constants in dic keys!
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
                '''
                print ':'.join([chromo, juc_start, juc_end])
                print juc_coord_dic[':'.join([chromo, juc_start, juc_end])]
                print(strand)
                print(i, j)
                print((juc_start, juc_end))
                raw_input(coords)
                '''
                juc_count += 1
    return juc_coord_dic

def index_annotations(annot_file, eventtype='SE'):
    '''
    Open miso annotation file and create indexed dictionaries.
    
    Dictionary format: 
        key: colon joined chromo, juc_start, juc_end
        value: dictionary with keys: eventid, start, end, chromo, strand, type 
        
    In this case, juc_start means END of an exon.
                  juc_end means START of the next exon.
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

def annotate_alexa_file(alexa_bed_file, junc_dic, output_file):
    '''
    With an indexed annotation of junctions read from MISO annotation, we will
    now read the alexa bed file, for each junction in alexa bed, we will
    search with the junction corresponds to any junctions from the indexed
    MISO annotation.
    
    If alexa junc matches with miso annotation, write a file containing
    the event details as well as the alexa junc id.
    
    Note, we assume junc_dic to have keys of form:
    "chromosome:intron_start:intron_end"
        We will rebuild this form and search the junc_dic if key exists.
    '''
    # Define colname indices in alexa bed file
    chromo_ind = 0
    alexa_start_ind = 1    # Start of junction, some bp upstrm from exon end.
    alexa_end_ind = 2    # End of junction, some bp dwnstrm of exon start.
    alexa_id_ind = 3
    blocksize_ind = 10
    # blockstarts_ind = 11
    
    # Define constant
    juc_str = 'juc'    # Alexa ID should start with this.
    
    # Define expected junc_dic subkey strings.
    eventid_str = 'eventid'
    start_str = 'start'
    end_str = 'end'
    chromo_str = 'chromo'
    strand_str = 'strand'
    type_str = 'type'
    
    # Define output header constants
    juc_id_str = 'junction_id'
    blocksizes_str = 'blocksizes'
    
    # Define count constants
    readcount = 0
    writecount = 0 
    
    with open(alexa_bed_file, 'rb') as readfile:
        reader = csv.reader(readfile, delimiter='\t')
        # Open write file similarly
        writefile = open(output_file, 'wb')
        writer = csv.writer(writefile, delimiter='\t')
        # Write header
        writer.writerow([chromo_str, start_str, end_str, juc_id_str, blocksizes_str,
                         eventid_str, strand_str, type_str])
        for row in reader:
            '''
            Bed file contains psr and juc, but psr only has 6 columns.
            Whereas juc has 11 columns. Use this fact to skip rows 
            that are not jucs.
            '''
            if len(row) <= 6:
                continue
            # Get data from row
            alexa_id = row[alexa_id_ind]
            # Check that alexa_id contains 'juc'
            if alexa_id[0:3] != juc_str:
                print('Warning, %s is not a junction alexa id.' %alexa_id)
            chromo = row[chromo_ind]
            alexa_start = row[alexa_start_ind]
            alexa_end = row[alexa_end_ind]
            # blocksize and starts are comma separated. So split'em.
            blocksize = row[blocksize_ind].split(',')
            # blockstarts = row[blockstarts_ind].split(',')
            '''
            Intron starts at end of exon.
            This means it can be calculated by 
            intron_start = alexa_start + size of first block.
            
            Intron ends at start of next exon.
            Calculate using:
            intron_end = alexa_start + 2nd element in block start
            or
            intron_end = alexa_end - size of 2nd block.
            '''
            try:
                intron_start = int(alexa_start) + int(blocksize[0])
                intron_end = int(alexa_end) - int(blocksize[1])
                # intron_end = int(alexa_start) + int(blockstarts[1])
            except ValueError:
                print('Could not convert one of the following into an integer:')
                for s in [alexa_start, alexa_end, blocksize[0], blocksize[1]]:
                    print s
            
            # Rebuild index and search dictionary for key.
            dic_key = ':'.join([chromo, str(intron_start), str(intron_end)])
            if dic_key in junc_dic:
                # Prepare list for writerow, initialize with
                # chromo, start, end, alexaid
                write_list = [chromo, alexa_start, alexa_end, alexa_id, 
                              row[blocksize_ind]]
                for subkey in [eventid_str, strand_str, type_str]:
                    '''
                    Values in subkey are a list, so join by comma.
                    We have iterated it so it goes eventid, strand, type.
                    This gives the row the proper order to match
                    the header.
                    '''
                    write_list.append(','.join(junc_dic[dic_key][subkey]))
                writer.writerow(write_list)
                writecount += 1
            else:
                pass
            readcount += 1
        writefile.close()
    return readcount, writecount

def main():
    if len(sys.argv) < 5:
        print('Some parameters need to be specified in command line.')
        sys.exit()
    alexa_bed_file = sys.argv[1]
    annot_file = sys.argv[2]    # From MISO
    output_file = sys.argv[3]
    eventtype = sys.argv[4]
    
    print('Indexing miso annotation for eventtype = %s...' %eventtype)
    junc_dic, _ = \
        index_annotations(annot_file, eventtype=eventtype)
    print('Created %s indexes.' %len(junc_dic.keys()))
    
    print(junc_dic.keys()[0:10])
    
    print('Reading alexaDB junctions and searching for AS events...')
    readcount, writecount = \
        annotate_alexa_file(alexa_bed_file, junc_dic, output_file)
    print('Lines read (# of alexa junctions): %s' %readcount)
    print('Lines written (# of ase events): %s' %writecount)
    
    
    

if __name__ == '__main__':
    main()
    