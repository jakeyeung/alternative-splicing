'''
Created on 2014-01-23

@author: jyeung
Parse constitutive exons, get flanking introns.
'''

import sys
import csv
from optparse import OptionParser

def get_chr_list(include_chr_prefix=False):
    '''
    Get list of chromosomes
    '''
    chr_list = [str(i) for i in range(1, 23)] + ['X', 'Y']
    if include_chr_prefix:
        chr_list = [''.join(['chr', i]) for i in chr_list]
    return chr_list

def get_start_end_strand(coord_id, include_chr_prefix=False):
    '''
    Coord id of form: chr7:95063412:95063861:-1
    Return chr, start, end ,strand
    that is: 7, 95063861, 95063412, -
    
    If remove_chr == True: returns 7
    if remove_chr == False: returns chr7
    '''
    coord_id_split = coord_id.split(':')
    chromo = coord_id_split[0]
    if include_chr_prefix:
        pass    # assumes default contains chr already
    elif not include_chr_prefix:
        # remove chr prefix from chromo
        chromo = chromo.split('chr')[-1]
    else:
        print 'include_chr_prefix must be True or False (boolean). \
            %s found.' %include_chr_prefix
    strand = coord_id_split[-1]
    # convert -1, +1 to - and +, respectively.
    if strand == '-1':
        strand = '-'
    elif strand == '+1':
        strand = '+'
    else:
        print 'Expected strand to be -1 or +1. %s found.' %strand
        sys.exit()
    # Get start and end by checking strand
    if strand == '+':
        start = coord_id_split[2]
        end = coord_id_split[1]
    else:
        start = coord_id_split[1]
        end = coord_id_split[2]
        
    # convert start, end to integer
    start = int(start)
    end = int(end)
    
    return chromo, start, end, strand

def get_flanking_introns(exon_start, exon_end, strand, length, location):
    '''
    Given exon start, end and strand, retrieve upstream and
    downstream flanking introns, set by user-defined length.
    
    Strand is important here.
    
    If location=='upstream', return:
        upstream_intron_start
        upstream_intron_end
    if location=='downstream', return:
        downstream_intron_start
        downstream_intron_end
    '''
    if strand == '+':
        if location=='upstream':
            # minus 1 to exclude first bp in exon start
            intron_end = exon_start - 1    
            intron_start = intron_end - length
        # add 1 to exclude last bp in exon_end
        elif location=='downstream':
            intron_start = exon_end + 1
            intron_end = intron_start + length
    elif strand == '-':
        if location=='upstream':
            # in negative strand, everything is reversed.
            intron_start = exon_start + 1
            intron_end = intron_start + length
        elif location=='downstream':
            intron_end = exon_end - 1
            intron_start = intron_end - length
    else:
        print 'Strand must be + or -, %s found.' %strand
        sys.exit()
    return intron_start, intron_end

def main():
    usage = 'usage: %prog input output_bed\n'\
        'Two args must be specified in commandline: \n'\
        '1) Input exons, ID in "chr7:95063412:95063861:-1" form.\n'\
        '2) Output bed file of flanking introns.\n'
    
    parser = OptionParser(usage=usage)
    
    parser.add_option('-c', '--column_index', dest='column_index',
                      default=0,
                      help='Column index containing ID, default 0')
    parser.add_option('-l', '--intron_length', dest='intron_length',
                      default=100,
                      help='How long of intron length to retrieve, '\
                        'default 100.')
    parser.add_option('-p', '--include_chr_prefix', dest='include_chr_prefix',
                      default=False,
                      help='Include chr prefix in chromosome name or not. '\
                        'True or False. Default False.')
    parser.add_option('--location', dest='location',
                      default='upstream',
                      help='Get upstream intron, downstream '\
                        'intron or exon from file.\n Must be upstream, '\
                            'downstream or exon. Default upstream')
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print 'Incorrect number of parameters specified.'
        print usage
        sys.exit()
    # parse args
    input_path = args[0]
    output_path = args[1]
    # BEGIN: parse options
    try:
        column_index = int(options.column_index)
    except ValueError:
        print '--column_index (-c) flag must be integer. %s found.' \
            %options.column_index
        sys.exit()
    try:
        length = int(options.intron_length)
    except ValueError:
        print '--intron_length (-l) flag must be integer. %s found.' \
            %options.intron_length
        sys.exit()
    include_chr_prefix = options.include_chr_prefix
    if include_chr_prefix in ['True', 'true', 'TRUE', 'T', True]:
        include_chr_prefix = True
    elif include_chr_prefix in ['False', 'false', 'FALSE', 'F', False]:
        include_chr_prefix = False
    else:
        print '--include_chr_prefix must be True or False. %s found.' \
            %include_chr_prefix
    location = options.location
    if location not in ['upstream', 'downstream', 'exon']:
        print '--location must be upstream, downstream or exon. %s found.' \
            %location
        sys.exit()
    # END: parse options
    
    print 'Output path: %s' %output_path
    print 'Location: %s' %location
    print 'Intron length (if applicable): %s' %length
    
    # init output
    output_file = open(output_path, 'wb')
    jwriter = csv.writer(output_file, delimiter='\t')
    
    # define bed score, a shading used for visualization
    bed_score = 999
    
    # define chromosomes, I want 1-22, X, Y only.
    chr_list = get_chr_list(include_chr_prefix=include_chr_prefix)
    
    # Read file
    writecount = 0
    with open(input_path, 'rb') as readfile:
        jreader = csv.reader(readfile, delimiter='\t')
        for rowcount, row in enumerate(jreader):
            '''
            coord_id of form: chr7:95063412:95063861:-1
            '''
            coord_id = row[column_index]
            # parse coord_id to get start, end
            chromo, exon_start, exon_end, strand = \
                get_start_end_strand(coord_id)
                
            # make sure chromosome is in chr_list
            if chromo not in chr_list:
                continue
            
            if location == 'exon':
                # start, ends are exon start ends
                start = exon_start
                end = exon_end
            elif location == 'upstream' or location == 'downstream':
                # get upstream and downstream introns from exon
                start, end= get_flanking_introns(exon_start, 
                                                 exon_end, 
                                                 strand, 
                                                 length, 
                                                 location)
            
            # create new ID to indicate upstream and downstream intron
            bed_id = ':'.join([coord_id, location])
            
            # write two rows per event, one upstream one downstream
            # write to file: chromosome, start, end, ID, score, strand
            row = [chromo, start, end, bed_id, bed_score, strand]
            jwriter.writerow(row)
            writecount += 1
    output_file.close()
    
    print '%s rows read, %s rows written to:\n%s' %(rowcount, writecount, output_path) 
    
if __name__ == '__main__':
    main()