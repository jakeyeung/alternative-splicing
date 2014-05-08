'''
Created on 2014-04-11

@author: jyeung
'''

from optparse import OptionParser
import sys

def main():
    usage = 'usage: %prog meme_gerp_summary_pkl\n'\
        'Requires one argument:\n'\
        '1) pkl file from summarize_meme_results'
    parser = OptionParser(usage=usage)
    (_, args) = parser.parse_args()
    
    if len(args) != 1:
        print usage
        sys.exit()
    
    fasta_path = args[0]
        
    in_frame = []
    out_frame = []
    in_frame_lengths = []
    out_frame_lengths = []
    seq_count = 0
    
    with open(fasta_path, 'rb') as fasta_file:
        for row in fasta_file:
            if row.startswith('>'):
                # get header ID
                fasta_header = row[1:]    # remove first character, '>'
            else:
                # check fasta_header is not None, ensures you have an ID
                if fasta_header is None:
                    print 'Error. No fasta header for row: %s' %row
                    sys.exit()
                # get sequence row, count length
                seq_length = len(row)
                if seq_length % 3 == 0:
                    # divisible by 3, save to "in-frame"
                    in_frame.append(fasta_header)
                    in_frame_lengths.append(seq_length)
                else:
                    # if not divisible by 3, save to out frame
                    out_frame.append(fasta_header)
                    out_frame_lengths.append(seq_length)
                # set fasta header to None
                fasta_header = None
                seq_count += 1
                
    # print how many are in frame, out frame.
    print in_frame_lengths
    print out_frame_lengths
    print '%s/%s in frame.' %(len(in_frame), seq_count)
    print '%s/%s out frame.' %(len(out_frame), seq_count)
    
    
if __name__ == '__main__':
    main()