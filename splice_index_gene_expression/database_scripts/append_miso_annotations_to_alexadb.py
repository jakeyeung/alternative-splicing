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

def create_juc_blocks_from_coords(coords, strand, eventtype='SE'):
    '''
    Strategy:
    1. Create "blocks" representing the two exons that are joined in a
        junction in upstream inclusion, exclusion and downstream inclusion (SE)
    2. Return "blocks" as a dictionary.
     
    '''
    # Init output dic
    juc_blocks_dic = {}
    # Init constant strings
    upstream_str = 'upstrm_juc_block'
    downstream_str = 'downstrm_juc_block'
    exclusion_str = 'exclusion_juc_block'
    keynames = [upstream_str, exclusion_str, downstream_str]
    
    coords = sorted(coords)    # Can differ depending on strand, so we sort.
    if eventtype=='SE':
        if strand=='+':
            '''
            We will define junctions by a list of tuples, each tuple is a 
            start-end of an exon body. Therefore, the list of tuples are
            the two exon bodies that the junction joins. 
            
            coords is sorted from smallest to largest, therefore, 
             - when strand is POSITIVE, upstream exon is (coords[0], coords[1])
             - when strand is NEGATIVE, upstream exon is (coords[4], coords[5])
              - coords should be a length of 6 because we are dealing with 
                  skipped exons.
            '''
            if len(coords) != 6:
                print('Length of coordinates list is not 6. List: %s' %coords)
                sys.exit()
            upstream_juc = [(coords[0], coords[1]), (coords[2], coords[3])]
            exclusion_juc = [(coords[0], coords[1]), (coords[4], coords[5])]
            downstream_juc = [(coords[2], coords[3]), (coords[4], coords[5])]
        elif strand=='-':
            downstream_juc = [(coords[0], coords[1]), (coords[2], coords[3])]
            exclusion_juc = [(coords[0], coords[1]), (coords[4], coords[5])]
            upstream_juc = [(coords[2], coords[3]), (coords[4], coords[5])]
        else:
            print('Expected %s strand to be + or -' %strand)
            sys.exit()
        # Store values into dictionary.
        juc_blocks_dic[upstream_str] = upstream_juc
        juc_blocks_dic[downstream_str] = downstream_juc
        juc_blocks_dic[exclusion_str] = exclusion_juc
        return juc_blocks_dic, keynames
    
def write_juc_blocks_to_dic(index_dic, eventid, chromo, exon_block_1,
                            exon_block_2, strand, juc_type):
    '''
    Given a junction that covers two exon blocks (block1 and block2), 
    create an indexed key and put relevant information such as 
    juc_type and strand into the dictionary. 
    
    Assumes exon_block_i is a tuple with:
        index 0 -> start
        index 1 -> end
    '''
    # Init subkeys
    eventid_str = 'eventid'
    block_1_start_str = 'block_1_start'
    block_1_end_str = 'block_1_end'
    block_2_start_str = 'block_2_start'
    block_2_end_str = 'block_2_end'
    chromo_str = 'chromo'
    strand_str = 'strand'
    type_str = 'type'
    block_start_str = 'block_start'
    block_end_str = 'block_end'
    
    key_block_start = ':'.join([chromo, exon_block_1[0], block_start_str])
    key_block_end = ':'.join([chromo, exon_block_1[1], block_end_str])
    
    for dic_key in [key_block_start, key_block_end]:
        if dic_key not in index_dic:
            # Initialize values as list so it is extendable. 
            index_dic[dic_key] = {eventid_str: [eventid],
                                  block_1_start_str: [int(exon_block_1[0])],
                                  block_1_end_str: [int(exon_block_1[1])],
                                  block_2_start_str: [int(exon_block_2[0])],
                                  block_2_end_str: [int(exon_block_2[1])],
                                  chromo_str: [chromo],
                                  strand_str: [strand],
                                  type_str: [juc_type]}
        else:
            '''
            If key already exists, append to list.
            '''
            for subkey, subval in zip\
                ([eventid_str, block_1_start_str, block_1_end_str, 
                  block_2_start_str, block_2_end_str, 
                  chromo_str, strand_str, type_str],
                 [eventid, int(exon_block_1[0]), int(exon_block_1[1]), 
                  int(exon_block_2[0]), int(exon_block_2[1]), 
                  chromo, strand, juc_type]):
                index_dic[dic_key][subkey].append(subval)
    return index_dic

def join_jucs_from_coords(coords, strand, eventtype='SE'):
    if eventtype=='SE':
        
        '''
        Strategy
        1. Order coordinates from smallest to largest (regardless of strand)
        2. Create three junction coordinates (upstream inclusion, 
            downstream inclusion, exclusion)
            2a. Upstream inclusion junction is:
                 - second (1) to third (2) coordinate for +
                 - fourth (3) to fifth (4) coordinate for -
            2b. Downstream inclusion junction is:
                 - fourth (3) to fifth (4) coordinate for +
                 - second (1) to third (2) coordinate for -
            2c. Exclusion junction is:
                 - second (1) to fifth (4) coordinate for + and -.
        3. Return junction starts and ends as a list in the following order:
            upstream (start, end), exclusion (start, end), downstream
        '''
        # Define keyname constants
        # exclusion_key = 'exclusion'
        # inclusion_up_key = 'inclusion_upstream'
        # inclusion_down_key = 'inclusion_downstream'
        # coord_dic = {}
        # Check there are 6 coordinates
        if len(coords) != 6:
            print('Expected 6 coordinates for an SE event.')
            sys.exit()
        # 1.
        coords = sorted(coords)
        # Get exclusion junction, strand independent.
        excl_juc_start = coords[1]
        excl_juc_end = coords[4]
        # Get inclusion junctions, which is strand dependent.
        if strand=='+':
            up_incl_juc_start = coords[1]
            up_incl_juc_end = coords[2]
            dwn_incl_juc_start = coords[3]
            dwn_incl_juc_end = coords[4]
        elif strand=='-':
            dwn_incl_juc_start = coords[1]
            dwn_incl_juc_end = coords[2]
            up_incl_juc_start = coords[3]
            up_incl_juc_end = coords[4]
        else:
            print('Expected strand to be + or -, %s found.' %strand)
            sys.exit()
        # Create list to summarize results.
        juc_coords_list = []
        for start, end in zip([up_incl_juc_start, 
                               excl_juc_start, 
                               dwn_incl_juc_start], 
                              [up_incl_juc_end, 
                               excl_juc_end, 
                               dwn_incl_juc_end]):
            juc_coords_list.append((start, end)) 
        return juc_coords_list

def get_coords_from_eventid(eventid, chromo, strand, 
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
    return coords, eventid_parsed_id
         
    '''
    Old crappy code...
      
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
                #print(strand)
                #print(i, j)
                #print((juc_start, juc_end))
                #raw_input(coords)
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
                #print ':'.join([chromo, juc_start, juc_end])
                #print juc_coord_dic[':'.join([chromo, juc_start, juc_end])]
                #print(strand)
                #print(i, j)
                #print((juc_start, juc_end))
                #raw_input(coords)
                juc_count += 1
    return juc_coord_dic
        '''

def index_annotations(annot_file, eventtype='SE'):
    '''
    Open miso annotation file (GFF3) and create indexed dictionaries.
    
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
    
    # String constants
    upstrm_incl_str = 'inclusion_upstream'
    excl_str = 'exclusion'
    dwnstrm_incl_str = 'inclusion_downstream'
    
    # Initialize dics
    index_dic = {}
    
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
        #TODO checking SE does nothing at the moment. 
        pass
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
            else:    # if == gene
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
                coords, eventid = \
                    get_coords_from_eventid(eventid, chromo, strand)
                juc_blocks_dic, keynames = \
                    create_juc_blocks_from_coords(coords, strand, eventtype)
                # Keynames are in order upstream, exclusion, downstream. 
                for key, juc_type in zip(keynames, ['inclusion_upstream', 
                                                    'exclusion', 
                                                    'inclusion_downstream']):
                    exon_block_1 = juc_blocks_dic[key][0]    # Tuple
                    exon_block_2 = juc_blocks_dic[key][1]
                    index_dic = write_juc_blocks_to_dic(index_dic, eventid, 
                                                        chromo, exon_block_1, 
                                                        exon_block_2, strand, 
                                                        juc_type)
                rowcount += 1
    return index_dic, rowcount

def annotate_alexa_file(alexa_bed_file, index_dic, output_file):
    '''
    With an indexed annotation of junctions read from MISO annotation, we will
    now read the alexa bed file, for each junction in alexa bed, we will
    search with the junction corresponds to any junctions from the indexed
    MISO annotation.
    
    If alexa junc matches with miso annotation, write a file containing
    the event details as well as the alexa junc id.
    
    Note, we assume index_dic to have keys of two forms:
    "chromosome:exon_start_pos:block_start" or 
    "chromosome:exon_end_pos:block-end"
        We will rebuild this form and search the index_dic if either key exists.
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
    
    # Define expected index_dic subkey strings.
    eventid_str = 'eventid'
    start_str = 'start'
    end_str = 'end'
    chromo_str = 'chromo'
    strand_str = 'strand'
    type_str = 'type'
    mismatch_score_str = 'mismatch_score'
    
    block_start_str = 'block_start'    # For recreating dic key
    block_end_str = 'block_end'    # For recreating dic key
    block_2_start_str = 'block_2_start'    # For searching subkey
    block_2_end_str = 'block_2_end'    # For searching subkey
    
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
            # Check that alexa_id contains 'juc' in first three letters.
            if alexa_id[0:3] != juc_str:
                print('Warning, %s is not a junction alexa id.' %alexa_id)
            chromo = row[chromo_ind]
            alexa_start = row[alexa_start_ind]
            alexa_end = row[alexa_end_ind]
            # blocksize and starts are comma separated. So split'em.
            blocksize = row[blocksize_ind].split(',')
            # blockstarts = row[blockstarts_ind].split(',')
            '''
            Get start and ends of first block/exon.
            block_1_start = alexa_start
            block_1_end = alexa_start + size of 1st block. 
            '''
            try:
                block_1_start = int(alexa_start)
                block_1_end = int(alexa_start) + int(blocksize[0])
                block_2_start = int(alexa_end) - int(blocksize[1])
                block_2_end = int(alexa_end)
                # block_1_end = int(alexa_start) + int(blockstarts[1])
            except ValueError:
                print('Could not convert one of the following into an integer:')
                for s in [alexa_start, alexa_end, blocksize[0], blocksize[1]]:
                    print s
            
            # Rebuild index and search dictionary for key.
            # Build a start key and an end key...
            block_start_key = ':'.join([chromo, str(block_1_start), 
                                        block_start_str])
            block_end_key = ':'.join([chromo, str(block_1_end), block_end_str])
            for dic_key in [block_start_key, block_end_key]:
                if dic_key in index_dic:
                    '''
                    # OK so either start or end matched, but make sure it also
                    # matches EITHER block_2_start OR block_2_end.
                    
                    If we can find the index/indices at which the block_2_start
                    or block_2_end matches, then we can write all the other 
                    information such as eventid, block start/ends 
                    and junction type into file. 
                    '''
                    index_list = []    # Initialize
                    for subkey, block_2_pos in zip([block_2_start_str, 
                                                    block_2_end_str], 
                                                    [block_2_start, 
                                                    block_2_end]):
                        # Intermediate list since we iterate twice.
                        temp_index_list = \
                        ([i for i, x in enumerate(index_dic[dic_key][subkey]) \
                            if x == block_2_pos])    
                        for temp_i in temp_index_list:
                            index_list.append(temp_i)
                    # Remove duplicates.
                    index_list = list(set(index_list))
                    
                    # If length not zero, then we write info to file.
                    if len(index_list) != 0:                    
                        # Prepare list for writerow, initialize with
                        # chromo, start, end, alexaid and blocksizes
                        write_list = [chromo, alexa_start, alexa_end, alexa_id, 
                                      row[blocksize_ind]]
                        # Add the MISO annotations now.
                        for subkey in [eventid_str, strand_str, type_str]:
                            '''
                            Values in subkey are a list, so join by comma.
                            We have iterated it so it goes eventid, strand, type.
                            This gives the row the proper order to match
                            the header.
                            '''
                            val_list = []
                            for i in index_list:
                                val_list.append(index_dic[dic_key][subkey][i])
                            write_list.append(','.join(val_list))
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
    
    print('Indexing miso annotation for eventtype == %s...' %eventtype)
    junc_dic, _ = \
        index_annotations(annot_file, eventtype=eventtype)
    print('Created %s indexes.' %len(junc_dic.keys()))
    
    # print(junc_dic.keys()[0:10])
    for k in junc_dic.keys()[0:20]:
        print k, junc_dic[k]
    
    print('Reading alexaDB junctions and searching for AS events...')
    readcount, writecount = \
        annotate_alexa_file(alexa_bed_file, junc_dic, output_file)
    print('Lines read (# of alexa junctions): %s' %readcount)
    print('Lines written (# of ase events): %s' %writecount)
    
if __name__ == '__main__':
    main()
    