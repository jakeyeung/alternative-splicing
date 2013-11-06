'''
Created on 2013-07-25

@author: jyeung
From MISO summary or differentially expressed isoforms, create a bed file
to be able to view the results on the genome browser. 
'''


import csv
import sys


def extract_chr_locs_from_summary(input_file, bed_output_file, trackname, 
                                  trackdescription):
    '''
    Read the summary file
    Read header
    Find chrom, strand, diff, start, stops in each row.
    Use abs(differential_psi) to determine the shading (higher shading meaning 
    higher abs(diff)
    Write to .bed file.
    '''
    print('Reading file %s...' %input_file)
    # Define column name constants
    chrom_str = 'chrom'
    strand_str = 'strand'
    start_str = 'mRNA_starts'
    end_str = 'mRNA_ends'
    diff_str = 'diff'
    event_str = 'event_name'
    header_str = 'track name=%s description="%s" useScore=1\n' %(trackname, 
                                                                 trackdescription)
    # Open files (didnt use with because python2.6 doesnt 
    # support two files using with)
    readfile = open(input_file, 'rb')
    writefile = open(bed_output_file, 'wb')
    # Write header
    writefile.write(header_str)
    # Create reader and writer objects.
    reader = csv.reader(readfile, delimiter='\t')
    writer = csv.writer(writefile, delimiter='\t')
    rheader = reader.next()
    
    rowcount = 0
    for row in reader:
        # Get chrom and strand
        chrom = row[rheader.index(chrom_str)]
        strand = row[rheader.index(strand_str)]
        diff = row[rheader.index(diff_str)]
        event = row[rheader.index(event_str)]
        start_full = row[rheader.index(start_str)]
        end_full = row[rheader.index(end_str)]
        
        # Convert diff to score
        score = abs(float(diff)) * 1000    # Convert scale between 0 to 1000.
        # Get event name, concatenate with diff.
        event_diff = ''.join([event, '_psi_diff:', diff])
        # Start and end is actually like a tuple (start1, start2), (end1, end2)
        # we will need to parse it out to get the start and end.
        # in cassette exons, start1 = start2, end1 = end2. 
        # Split string by commas, giving a list of ['(start1', 'start2)']
        # Take first in list, remove the beginning bracket, then convert to int.
        start = int(start_full.split(',')[0][1:])
        end = int(end_full.split(',')[0][1:])
        # Write info to bed file.
        writer.writerow([chrom, start, end, event_diff, score, strand])
        rowcount += 1
    readfile.close()
    writefile.close()
    print('%s rows read and written to %s.' %(rowcount, bed_output_file))

def main():
    if len(sys.argv) < 3:
        print('Summary file and output file must be specified in command line.')
        sys.exit()
    summary_file = sys.argv[1]
    bed_output_file = sys.argv[2]
    trackname = sys.argv[3]    # no spaces
    trackdescription = sys.argv[4]    # underscore separated words.
    
    # break track description into separate words.
    trackdescription = trackdescription.replace('_', ' ')
    
    extract_chr_locs_from_summary(summary_file, bed_output_file, trackname, 
                                  trackdescription)

if __name__ == '__main__':
    main()