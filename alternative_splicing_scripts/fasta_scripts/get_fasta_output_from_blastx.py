'''
Created on 2013-11-22

@author: jyeung
From blastx output (from command line and option 
flag -outfmt "6 qseqid qlen length pident qseq"),
Extract the non-redundant 100% matches and get their
amino acid sequence.
Output the 100% matches to a fasta file.
'''

import sys
import csv
from optparse import OptionParser

def index_blastx_output(blastx_output_path, qseqid_index, 
                                     pident_index, qseq_index):
    '''
    Stores blastx output as a dictionary.
    Dictionary contains:
    {qseqid: qseq}
    Only stores values for pident = 100% and 
    no redundancies. We will pick the first 100% and ignore
    others that match.
    '''
    blastx_dic = {}
    with open(blastx_output_path, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        # Expects no headers, so iterate rows right away.
        for row in myreader:
            qseqid = row[qseqid_index]
            pident = float(row[pident_index])
            qseq = row[qseq_index]
            if pident == 100:
                blastx_dic[qseqid] = qseq
    return blastx_dic

def write_fasta_sequence(seqid, seq, writefile, max_length=70):
    '''
    Write to writefile (not csv obj) the seqid and seq in
    fasta format.
    Example:
    >seqid
    VFRVYQGQQPGTCM
    Make new line every 70 characters...
    '''
    # Write description
    writefile.write('>%s\n' %seqid)
    # init counters
    letter_count = 0
    current_seq = []
    for s in seq:
        current_seq.append(s)
        letter_count += 1
        if letter_count % max_length == 0:
            # 70 letters collected, write as a line.
            writefile.write(''.join(current_seq + ['\n']))
            # Reset my counters
            current_seq = []
    # Finished looping sequence, write remaining sequences
    # but first check current_seq is not empty.
    if len(current_seq) != 0:
        writefile.write(''.join(current_seq + ['\n']))
    return None
            
def write_blastx_dic_to_file(blastx_dic, fasta_path):
    '''
    Makes a fasta file with seqid as id
    Sequences will be 70 characters long before forming new line.
    '''
    # Init fasta writer
    writecount = 0
    with open(fasta_path, 'wb') as writefile:
        for seqid in blastx_dic.keys():
            seq = blastx_dic[seqid]
            write_fasta_sequence(seqid, seq, writefile)
            writecount += 1
    print '%s sequences written to: %s' %(writecount, fasta_path)
    return writecount

def main():
    usage = 'usage: %prog [options] blastx_output.out fasta_file.fasta \n'\
        'Two args must be specified in commandline: \n'\
        '1) output from blastx, in tabular format.\n'\
        '2) fasta file output (to be written).\n'
    parser = OptionParser(usage=usage)
    parser.add_option('-i', '--qseqid_index', dest='qseqid_index',
                      help='Column index (0-based) containing sequence ID '\
                      '\nExample, if blastx -outfmt option was'\
                      '"6 qseqid qlen length pident qseq", -i should be 0 '\
                      '(1st column). Default is 0',
                      default=0)
    parser.add_option('-p', '--pident_index', dest='pident_index',
                      help='Column index (0-based) containing percentage of '\
                      'identical matches.\nExample, if blastx -outfmt option was'\
                      '"6 qseqid qlen length pident qseq", -p should be 3 '\
                      '(4th column). Default is 3',
                      default=3)
    parser.add_option('-s', '--qseq_index', dest='qseq_index',
                      help='Column index (0-based) amino acid sequence '\
                      '.\nExample, if blastx -outfmt option was'\
                      '"6 qseqid qlen length pident qseq", -p should be 4 '\
                      '(5th column). Default is 4',
                      default=4)
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print(usage)
        sys.exit()
    blastx_output_path = args[0]
    fasta_file_path = args[1]
    qseqid_index = options.qseqid_index
    pident_index = options.pident_index
    qseq_index = options.qseq_index
    
    blastx_dic = index_blastx_output(blastx_output_path, qseqid_index, 
                                     pident_index, qseq_index)
    write_blastx_dic_to_file(blastx_dic, fasta_file_path)

if __name__ == '__main__':
    main()