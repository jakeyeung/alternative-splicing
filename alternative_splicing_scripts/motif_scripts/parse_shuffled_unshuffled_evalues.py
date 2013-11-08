'''
Created on 2013-11-06

@author: jyeung
After running meme_batch_mode.sh for shuffled and unshuffled sequences,
we want to get the evalues from both shuffled and unshuffled sequences.

Parse through all the textfiles in shuffled and unshuffled. Then write to textfile
that is then used for plotting in R.

Assumes meme output filename is meme.txt
'''

import sys
import os
import csv
from optparse import OptionParser
from utilities import get_evalues

def main():
    usage = 'usage: %prog motif_dir shuffled_dir output_file.parsed'
    parser = OptionParser(usage=usage)
    (_, args) = parser.parse_args()
    if len(args) < 3:
        print('There args must be specified in commandline: \nMotif Directory '\
              '(unshuffled) \nMotif Directory (shuffled) \nOutput file')
        sys.exit()
    
    non_null_rootdir = args[0]
    null_rootdir = args[1]
    write_out_file = args[2]
    
    # Get list of directories containing meme outputs.
    non_null_dirs = os.listdir(non_null_rootdir)
    null_dirs = os.listdir(null_rootdir)
    
    # Init write out file
    writefile = open(write_out_file, 'wb')
    mywriter = csv.writer(writefile, delimiter='\t')
    
    # Loop non_null_files and null_files, grabbing their evalues.
    non_null_evalues = []
    null_evalues = []
    for non_null_dir, null_dir in zip(non_null_dirs, null_dirs):
        # Get meme file from directory. Assumes filename as meme.txt
        non_null_file = os.path.join(non_null_rootdir, 
                                     non_null_dir, 
                                     'meme.txt')
        null_file = os.path.join(null_rootdir, 
                                 null_dir, 
                                 'meme.txt')
        non_null_evalues += get_evalues.retrieve_evalues_from_file(non_null_file)
        null_evalues += get_evalues.retrieve_evalues_from_file(null_file)
    
    # Write evalues to file, into a format that is easy to ggplot2 in R.
    # write header
    mywriter.writerow(['Group', 'Motif E-Value'])
    writecount = 0
    for non_null_eval, null_eval in zip(non_null_evalues, null_evalues):
        # Call non-nulls "Sequences of Interest"
        # Call nulls "Shuffled Sequences"
        mywriter.writerow(['Sequences of Interest', non_null_eval])
        mywriter.writerow(['Shuffled Sequences', null_eval])
        writecount += 2
    
    print '%s evalues written to file: %s' %(writecount, write_out_file)
    writefile.close()

if __name__ == '__main__':
    main()