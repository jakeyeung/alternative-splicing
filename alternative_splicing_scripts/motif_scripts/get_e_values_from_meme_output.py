'''
Created on 2013-10-26

@author: jyeung

Parses meme output and retrieves evalues for each of the 10 motifs.
This is currently written to work for MEME outputs via webserver.
#TODO This should be optimized for running in commandline in the future.
'''

from optparse import OptionParser
import sys
import os
import csv
import re

def retrieve_evalues_from_file(meme_file):
    matchcount = 0
    evalues = []
    with open(meme_file, 'rb') as readfile:
        for line in readfile:
            m = re.search('(?<=E-value = )\S+', line)
            if m:
                matchcount += 1
                evalues.append(m.group(0))
    if matchcount != 10:
        print 'Warning: expected 10 (%s found) Evalues to be retrieved from %s'\
             %(matchcount, meme_file)
    return evalues

def get_full_path_list(mydir, ext=None):
    '''
    Given directory, retrieve all files with extension.
    Return file as a list of full paths.
    '''
    myfiles = [f for f in os.listdir(mydir) if f.endswith(ext)]
    myfiles_fullpaths = [os.path.join(mydir, f) for f in myfiles]
    return myfiles_fullpaths

def main():
    usage = 'usage: %prog non_null_directory null_directory output_file.parsed'
    parser = OptionParser(usage=usage)
    (_, args) = parser.parse_args()
    if len(args) < 3:
        print('Two directories (one non-null, one null) and output file must be'\
              ' specified in commandline.')
        sys.exit()
    
    non_null_dir = args[0]
    null_dir = args[1]
    write_out_file = args[2]
    
    # Get list of .out files from directory.
    '''
    non_null_files = \
        [f for f in os.listdir(non_null_dir) if f.endswith('.out')]
    null_files = \
        [f for f in os.listdir(null_dir) if f.endswith('.out')]
    '''
    non_null_files = get_full_path_list(non_null_dir, ext='out')
    null_files = get_full_path_list(null_dir, ext='out')
        
    # Init write out file
    writefile = open(write_out_file, 'wb')
    mywriter = csv.writer(writefile, delimiter='\t')
    
    # Loop non_null_files and null_files, grabbing their evalues.
    non_null_evalues = []
    null_evalues = []
    for f, null_f in zip(non_null_files, null_files):
        non_null_evalues += retrieve_evalues_from_file(f)
        null_evalues += retrieve_evalues_from_file(null_f)
    
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