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
from append_multiple_bed_files import append_multiple_bed_files

def get_chrom_from_event(event_name):
    '''
    From event name, parse the string to get chromosome.
    Example:
    chr9:131772050:131772150:-@chr9:131771732:131771746:-@chr9:131771385:131771625:-
    Get chr9 (first element in ':' split string)
    '''
    chromo = event_name.split(':')[0]
    return chromo

def get_strand_from_event(event_name):
    '''
    Parse event name string to get strand.
    Example:
    chr9:131772050:131772150:-@chr9:131771732:131771746:-@chr9:131771385:131771625:-
    Get - (last element in ':' split string)
    '''
    strand = event_name.split(':')[-1]
    return strand 

def create_exon_coords(event_name):
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
    the second (start) and third (end) element of each group.
    Ignore first element because it is chrom. 
    
    NOTE: UCSC uses 0-base start and 1-based end! Therefore,
    we have to subtract ONE from the exon-start.
    '''
    event_name_list = event_name.split(':')
    # Loop through event_name_list, get every 2nd and 3rd element
    # as a start/end pair.
    exon_starts = []
    exon_ends = []
    for i, coordinate, in enumerate(event_name_list):
        if i % 3 == 0:
            pass    # no remainder means you are at first element/chrom
        elif i % 3 == 1:
            exon_starts.append(int(coordinate)-1)    # 
        elif i % 3 == 2:
            exon_ends.append(int(coordinate))
        else:
            sys.exit('Parsing error in %, %' %(i, coordinate))
    
    return exon_starts, exon_ends

def create_exon_intron_bed_files(output_bed_file, bed_description, eventtype='SE'):
    '''
    Create dictionary containing writefile objects depending on 
    event type.
    
    This allows writing of many files using loops. 
    
    Also writes bed_header to file. 
    '''
    if eventtype == 'SE':
        # Initialize list of indices
        exon_index = range(1, 4)    # 3 exons (1, 2, 3)
        intron_index = range(1, 3)    # 2 introns (1, 2)
        
        # Create keynames for exon and intron dictionaries.
        exon_keys = [''.join(['exon_', str(i)]) for i in exon_index]
        intron_keys = [''.join(['intron_', str(i)]) for i in intron_index]
        
        # Record the full path names used to be returned.
        fullpaths = []
        
        # Create dic to store output write objects to write exon coords
        output_file_dic_exons = {}
        for i, k in zip(exon_index, exon_keys):
            # Append output_bed_file suffix
            output_bed_file_suffixed = \
                ''.join([output_bed_file, '_exon', str(i), '.bed'])
            fullpaths.append(output_bed_file_suffixed)
            output_file_dic_exons[k] = open(output_bed_file_suffixed, 'wb')
        
        # Create dic to store write objects to write intron coords
        output_file_dic_introns = {}
        for i, k in zip(intron_index, intron_keys):
            # Append output_bed_file suffix
            output_bed_file_suffixed = \
                ''.join([output_bed_file, '_intron', str(i), '.bed'])
            fullpaths.append(output_bed_file_suffixed)
            output_file_dic_introns[k] = open(output_bed_file_suffixed, 'wb')
            
        # Writes bed header to file
        for k in exon_keys:
            write_str = 'track name=%s description="%s"\n' %(k, bed_description)
            output_file_dic_exons[k].write(write_str)
        for k in intron_keys:
            write_str = 'track name=%s description="%s"\n' %(k, bed_description)
            output_file_dic_introns[k].write(write_str)
        
    else:
        print('Unknown event type (or not yet implemented.')
        sys.exit()
        
    return output_file_dic_exons, output_file_dic_introns, fullpaths

def create_intron_coords_from_exon_coords(exon_starts, exon_ends, strand, eventtype='SE'):
    '''
    Given exon_starts, exon_ends, create its corresponding intron_start
    and intron_ends.
    
    Input:
        exon_starts: list of coordinates of the beginning of each exon.
        exon_ends: list of coordinates of the END of each corresponding exon.
        strand: positive or negative (important to determine whether
        exon starts are increasing [+] or decreasing [-]. 
        
        Example:
        First element in start is start of first exon.
        First element in end is end of first exon.
        
    Output:
        intron_starts: list of coordinates of the beginning of intron (end of exon).
        intron_ends: list of coordinates of the END of the intron (beginning
        of next exon).
    ''' 
    '''
    Method for positive strand:
    # Introns will not include the first exon_start nor the last exon_end, because
    # introns only connect exons. First and last exon have nothing to connect.
    # 
    # Therefore, we remove first element in exon_start and last element in exon_end,
    # then we connect the corresponding exons together and call them introns.
    
    Also:
        -intron starts at one bp AFTER end of exon.
        -intron ends at one bp BEFORE start of next exon.
    
    # NOTE:
    UCSC uses the zero-based start and the 1-based end!
    Therefore, intron-start was de-incremented by 1, 
        intron-end was not modified.
    '''
    '''
    Method for negative strand:
    exon start/ends are decreasing in order. 
    Intron 1 will therefore be
    defined as (skipped exon end + 1) to (exon1 start - 1)
    Intron 2 will be defined as (exon3 end + 1) to (skippedexon start - 1)
    '''
    # Define constants
    if eventtype == 'SE':
        # Take 300 bps around skipped exon instead of full intron.
        dist_around_se = 300    
    
    # Clip ends, but differently depending on strand + or strand -
    if strand == '+':
        exon_ends_clipped = exon_ends[:-1]    # Last exon end is 3' end
        exon_starts_clipped = exon_starts[1:]    # First exon end is 5' start
    elif strand == '-':
        exon_ends_clipped = exon_ends[1:]    # First exon_end is 5' start.
        exon_starts_clipped = exon_starts[:-1]    # Last exon_start is 3' end.
    # Get intron starts and ends  
    try:
        intron_starts = [int(i) for i in exon_ends_clipped]
        intron_ends = [int(i) for i in exon_starts_clipped]
    except ValueError:
        print('Could not increment integers in one or both of the lists:')
        print(exon_ends_clipped)
        print(exon_starts_clipped)
        sys.exit()
    
    '''    
    # For skipped exon events (SE), we will extract 300 bp around the 
    # skipped exon for the intron area.
    
    Use second exon (index 1) and get upstream and downstream 300 bp.
    
    Use that value as new first intron start or second intron end 
    only if it does not overlap with the upstream/downstream exon.
    '''
    if eventtype == 'SE':
        # Get upstream and downstream 300 bp around SE.
        upstrm_se = exon_starts[1] - dist_around_se
        dwnstrm_se = exon_ends[1] + dist_around_se
        if strand == '+':
            intron_ends[1] = min(dwnstrm_se, intron_ends[1])
            intron_starts[0] = max(upstrm_se, intron_starts[0])
        elif strand == '-':
            intron_ends[0] = min(dwnstrm_se, intron_ends[0])
            intron_starts[1] = max(upstrm_se, intron_starts[1])
            
    return intron_starts, intron_ends

def write_coords_to_multi_files(output_bed_file, output_file_dic_exons, 
                                output_file_dic_introns,
                                chrom, exon_starts, exon_ends, 
                                intron_starts, intron_ends, event,
                                score, strand):
    '''
    List of exon-starts and ends, write to multiple files, one for each
    exon and intron region. 
    '''
    # Check lengths
    if len(exon_starts) != len(exon_ends):
        print('Warning, length of exons starts and ends are not equal!'\
         '(%s vs. %s)' %(len(exon_starts), len(exon_ends)))
        print('exon_starts: %s' %exon_starts)
        print('exon_ends: %s' %exon_ends)
        raw_input('Continue? (Press enter:')
    
    '''
    # Create corresponding intron start and ends.
    intron_starts, intron_ends = \
        create_intron_coords_from_exon_coords(exon_starts, exon_ends)
    '''
    
    # Get the sorted exon_key and intron_keys for looping later...
    exon_keys = sorted(output_file_dic_exons.keys())
    intron_keys = sorted(output_file_dic_introns.keys())
    
    '''
    # Initialize list of indices
    exon_index = range(0, len(exon_starts))
    intron_index = range(0, len(exon_starts) - 1)
    
    # Create keynames for exon and intron dictionaries.
    exon_keys = [''.join(['exon_', str(i)]) for i in exon_index]
    intron_keys = [''.join(['intron_', str(i)]) for i in intron_index]
    
    # Create dic to store output write objects to write exon coords
    output_file_dic_exons = {}
    for i, k in zip(exon_index, exon_keys):
        # Append output_bed_file suffix
        output_bed_file_suffixed = \
            ''.join([output_bed_file, '_exon', str(i), '.bed'])
        output_file_dic_exons[k] = open(output_bed_file_suffixed, 'wb')
    
    # Create dic to store write objects to write intron coords
    output_file_dic_introns = {}
    for i, k in zip(intron_index, intron_keys):
        # Append output_bed_file suffix
        output_bed_file_suffixed = \
            ''.join([output_bed_file, '_intron', str(i), '.bed'])
        output_file_dic_introns[k] = open(output_bed_file_suffixed, 'wb')
    '''
    
    writecount = 0
    # Write exons to exon files.
    for start, end, k in zip(exon_starts, exon_ends, exon_keys):
        # Create tab separated string for writing to file.
        write_str = '\t'.join([chrom, str(start), str(end), event, 
                               score, strand, '\n'])
        output_file_dic_exons[k].write(write_str)
        writecount += 1
        
    # Write introns to intron files.
    for i_start, i_end, i_k in zip(intron_starts, intron_ends, 
                                   intron_keys):
        '''
        # introns connect an end to a start. 
        # We will connect first exon-end to second exon-start, all the way
        # until the connection of second last exon-end to last exon-start.
        # 
        # Therefore, we use all coordinates in exon_ends/starts EXCEPT:
            -exon_ends[0]
            -exon_starst[-1]
        '''
        write_str = '\t'.join([chrom, str(i_start), str(i_end), event, 
                               score, strand, '\n'])
        output_file_dic_introns[i_k].write(write_str)
        writecount += 1
        
    return writecount

def extract_coordinates_from_miso_bf(miso_file, output_bed_file, 
                                     bed_description, eventtype='SE'):
    '''
    Input:
        miso_file - file from bayes analysis. Expects column names such as
            event_name, chrom, and strand.
    Output:
        output_bed_file - bed file extracting coordinates of upstream, skipped
        and downstream exon as well as the two upstream/downstream introns.
        
        Each region (upstream skipped, downstream, up/down introns) are outputted
        into a separate bed file, with the output_bed_file suffixed to indicate
        the type of region.
    
    #TODO: chrom_str and strand_str can be obtained from event_name, so it 
    can be replaced by a function.
    '''
    # Define column name constants
    event_name_str = 'event_name'
    event_name_str2 = 'event'    # Another possible name.
    chrom_str = 'chrom'    # May not exist as colname
    strand_str = 'strand'    # May not exist as colname. 
    score = '999'    # For UCSC shading
    
    '''
    # TODO: delete me. 
    # Create read_write_obj for reading and writing simultaneously and also works
    # in python 2.6
    # read_write_obj = read_write(miso_file, output_bed_file, header=True)
    '''
    with open(miso_file, 'rb') as readfile:
        reader = csv.reader(readfile, delimiter='\t')
        # Create open file objects.
        if eventtype == 'SE':
            '''
            Create three exon bed files and two intron files,
            corresponding to:
                upstream exon
                skipped exon
                downstream exon
                upstream from skipped intron
                downstream from skipped intron
            '''
            output_file_dic_exons, output_file_dic_introns , full_outpaths = \
                create_exon_intron_bed_files(output_bed_file, 
                                             bed_description, 
                                             eventtype='SE')
        # Get header
        colnames = reader.next()
        wcount = 0
        for row in reader:
            try:
                event_name = row[colnames.index(event_name_str)]
            except ValueError:    # Try second string.
                event_name = row[colnames.index(event_name_str2)]
            try:
                chrom = row[colnames.index(chrom_str)]
            except ValueError:    # May not exist.
                chrom = get_chrom_from_event(event_name)
            try:
                strand = row[colnames.index(strand_str)]
            except ValueError:    # May not exist.
                strand = get_strand_from_event(event_name)
            
            # Get exon and intron coordinates
            exon_starts, exon_ends = \
                create_exon_coords(event_name)        
            intron_starts, intron_ends = \
                create_intron_coords_from_exon_coords(exon_starts, exon_ends, strand)
            
            writecount = \
                write_coords_to_multi_files(output_bed_file, 
                                            output_file_dic_exons, 
                                            output_file_dic_introns,
                                            chrom, exon_starts, exon_ends, 
                                            intron_starts, intron_ends,
                                            event_name,
                                            score, strand)
            wcount += writecount
        # Close exon and intron files.
        for writeobj in (output_file_dic_exons.values() + 
                         output_file_dic_introns.values()):
            writeobj.close()
             
        return wcount, full_outpaths

def main():
    if len(sys.argv) < 3:
        print('Miso_bf file (or bh_adj t-test summary file)'\
              ' and output_bed_file '\
              'must be specified in command line..')
        sys.exit()
    miso_file = sys.argv[1]
    output_bed_file = sys.argv[2]
    appended_bed_file_path = sys.argv[3]
    bed_description = sys.argv[4]
    
    '''
    If user inputted output_bed_file with a .* ending, we will remove it
    because we want to append more strings to end of filename before 
    tagging the .bed ending. 
    '''
    output_bed_file = '.'.join(output_bed_file.split('.')[:-1])
    wcount, outpaths = extract_coordinates_from_miso_bf(miso_file, 
                                                        output_bed_file, 
                                                        bed_description)
    append_multiple_bed_files(outpaths, appended_bed_file_path)
    print('Done. %s total coordinates extracted.' %wcount)
    print('Bed files written to:')
    for f in outpaths:
        print(f)
    print('Appended bed file written to:\n%s' %appended_bed_file_path)

if __name__ == '__main__':
    main()