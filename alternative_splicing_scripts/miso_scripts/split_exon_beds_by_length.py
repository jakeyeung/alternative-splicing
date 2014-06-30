'''
Created on 2013-12-02

@author: jyeung

Exon lengths are long. So let's split them by
a certain length that's user defined (e.g. 100 bp)

Care must be taken because we split bed files that
are both + and - strand.

We want 100 bp from start and end of exon, which differs
depending if it is + or - strand.
'''

import sys
import os
import csv
from optparse import OptionParser
from split_beds_into_inclusion_exclusion import modify_bed_header

def get_dic_keys_5p_3p():
    '''
    Generate dic keys for 5p and 3p.
    Order matters.
    '''
    return ['5p', '3p']
    
def add_5p_3p_to_filename(myfilename):
    '''
    Adds 5p and 3p to filename, returns 
    as a dictionary containing filepaths.
    '''
    # Define suffixes to add to output filename
    exon_dic_keys = get_dic_keys_5p_3p()    # ['5p', '3p']
    output_files_dic = {}
    
    for suffix in exon_dic_keys:
        # Strip extension, add suffix, then add back extension
        split_file = os.path.splitext(myfilename)
        output_files_dic[suffix] = ''.join([split_file[0], 
                                            '_', 
                                             suffix, 
                                             split_file[1]])
    return output_files_dic

def check_exon_lengths_assign_nbp(exon1, exon2, n=100):
    '''
    Check exon lengths, make sure they are not
    smaller than n. 
    If exon lengths > n, then add n to 5p, subtract n
    to 3p.
    
    If it is smaller, then you cannot add n because
    it would go into intron region. Rather, just take
    the opposite end because that's the farthest you can
    go without hitting intron region.
    
    n is default 100, meaning 100 bp upstream or downstream.
    
    If this is pos strand,
    exon1 is 5p.
    exon2 is 3p.
    '''
    splitcount = 0
    exon_len = abs(exon1 - exon2)
    if exon_len > n:
        exon1_nbp = exon1 + n
        exon2_nbp = exon2 - n
        splitcount += 1
    elif exon_len <= n:
        exon1_nbp = exon2
        exon2_nbp = exon1
    return exon1_nbp, exon2_nbp, splitcount

def iterate_rows_split_exons(jreader, jwriter_5p, jwriter_3p, 
                             length_to_split=100):
    '''
    # Iterate through rows, modify start and ends depending on strand,
    # then write new modified start-ends to respective 5p and 3p files.
    Bed file is of form:
    row[0]: chromosome
    row[1]: exon start if pos, exon end if neg
    row[2]: exon end if pos, exon start if pos
    row[3]: miso id
    row[4]: score
    row[5]: strand
    
    Length to split: default to 100.
    
    If exon length is shorter than 100, then do not adjust exon starts/ends.
    5p will equal 3p in that case.
    
    Readers and writers are all CSV write/read objects. 
    '''
    splitcount = 0
    for rowcount, row in enumerate(jreader):
        strand = row[5]
        # exon1 is exon start if pos, end if neg.
        # exon2 is exon end if pos, start if neg.
        exon1 = int(row[1])
        exon2 = int(row[2])
        exon1_nbp, exon2_nbp, scount = \
            check_exon_lengths_assign_nbp(exon1, 
                                          exon2, 
                                          length_to_split)
        splitcount += scount
        # Write new coordinates to respective files.
        if strand == '+':
            '''
            If positive strand:
                make row[2] exon1_nbp in 5p file
                make row[1] exon2_nbp in 3p file
            '''
            row_5p = row
            row_3p = row
            row_5p[2] = exon1_nbp
            row_3p[1] = exon2_nbp
            jwriter_5p.writerow(row_5p)
            jwriter_3p.writerow(row_3p)
        elif strand == '-':
            '''
            If negative strand:
                make row[1] exon2_nbp in 5p file
                make row[2] exon1_nbp in 3p file
            '''
            row_5p = row
            row_3p = row
            row_5p[1] = exon2_nbp
            row_3p[2] = exon1_nbp
            jwriter_5p.writerow(row_5p)
            jwriter_3p.writerow(row_3p)
    print '%s exons were split.' %splitcount
    return rowcount
            
def main():
    usage = 'usage: %prog [opts] exon_bed_file.bed\n'\
        'One arg must be specified in commandline: \n'\
        '1) Bed file containing full exons.\n--help for option help'
    parser = OptionParser(usage=usage)
    parser.add_option('-o', '--output_bed_name', 
                      dest='output_bed_name',
                      help='Output bed name,'\
                      'Default is exon_1_3pOR5p.bed (N is length of split)',
                      default='')
    parser.add_option('-l', '--length', dest='length',
                      help='Length of exon to split 5p and 3p. Default 100.',
                      default='100')
    (options, args) = parser.parse_args()
    
    if len(args) < 1:
        print usage
        sys.exit()
    exon_bed_file = args[0]
    
    # Define output filename inside dics
    # Filenames accessible by dic["5p"] or dic["3p"]
    if options.output_bed_name == '':
        # If default, append exon_bed_file to suffix as output.
        output_files_dic = add_5p_3p_to_filename(exon_bed_file)
    else:    # if not default, append user defined to suffix as output
        output_files_dic = add_5p_3p_to_filename(options.output_bed_name)
    
    # initialize write files
    writefile_5p = open(output_files_dic['5p'], 'wb')
    writefile_3p = open(output_files_dic['3p'], 'wb')
    
    # read bed file, for each exon (row), create 5p and 3p ends of each exon
    with open(exon_bed_file, 'rb') as readfile:
        # Read header, modify it then write to each writefile.
        readheader = readfile.readline()
        header_5p = modify_bed_header(readheader, '5p', append_start=10)
        header_3p = modify_bed_header(readheader, '3p', append_start=10)
        for writefile, writeheader in \
            zip([writefile_5p, writefile_3p], [header_5p, header_3p]):
            writefile.write(writeheader)
        # Create csv objects for read and writefiles, then iterate read rows
        myreader = csv.reader(readfile, delimiter='\t')
        mywriter_5p = csv.writer(writefile_5p, delimiter='\t')
        mywriter_3p = csv.writer(writefile_3p, delimiter='\t')
        rowcount = iterate_rows_split_exons(myreader, 
                                            mywriter_5p, 
                                            mywriter_3p, 
                                            options.length)
        for f in output_files_dic.values():
            print '%s rows written to file: %s' %(rowcount, f)
        
if __name__ == '__main__':
    main()