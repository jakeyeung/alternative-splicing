'''
Created on 2013-05-30

@author: jyeung
From the OUTPUT from script 'find_splice_location.py', set it up 
into an output text table to be used for OptDis. 

The readfile should be a table of either:
1) probe expressions and gene expressions
2) SI and gene expressions
With one gene, one data point per sample. 
'''


import os
import sys
from utilities import set_directories, read_write_gene_info


# Set constants

# Folder names
outputdir = 'output'
tabledir = 'tables'

# Input column names
gene_symbol_colname = 'gene_symbol'
# Sample constants
sample_list = ['X946_Urethra', 
                'X972.2.Penila',
                'AB352',
                'X963.Lpa.JP',
                'X963.L.LN',
                'X1005',
                'X890L',
                'X945',
                'X961']
classifier_dict = {'NEPC':['X946_Urethra', 'X972.2.Penila', 'AB352', 'X963.Lpa.JP', 'X963.L.LN',], 
                   'PC':['X1005', 'X890L', 'X945', 'X961']}
full_sample_length = 13    # All samples, including cell lines and lymph nodes.

# Output column names, only one defined here, the others are defined 
# inside the 'with' statement.
gene_label_colname = 'gene_name'

# Set directories
mydirs = set_directories.mydirs('', outputdir)

def extract_and_write_file(readwrite):
    # Read rows, extract gene information. 
    with readwrite:
        # Write column names
        outputcolnames = [gene_symbol_colname]
        for tumor_class, sample in classifier_dict.iteritems():
            for _ in sample:
                outputcolnames.append(tumor_class)
        readwrite.writenext(outputcolnames)
        '''
        Read and write probe or SI information. 
        '''
        while True:
            try:
                row = readwrite.readnext()
            except StopIteration:
                print('%s rows read, %s written. Done' %(readwrite.readrowcount, 
                                                         readwrite.writerowcount))
                break
            # Initialize writerow by putting gene symbol as first element.
            writerow = [row[readwrite.inputcolnames.index(gene_symbol_colname)]]
            # Now append probe or SI values to list.
            for s in sample_list:
                writerow.append(row[readwrite.inputcolnames.index(s)])
            # Write to file.
            readwrite.writenext(writerow)
    return None
                
if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Input and output filename must be provided in the command line.')
        sys.exit()
    read_fname = sys.argv[1]
    write_fname = sys.argv[2]
    
    read_path = os.path.join(mydirs.outputdir, tabledir, read_fname)
    write_path = os.path.join(mydirs.outputdir, tabledir, write_fname)
    
    readwrite = read_write_gene_info.read_write_gene_info(read_path, write_path)
    
    extract_and_write_file(readwrite)