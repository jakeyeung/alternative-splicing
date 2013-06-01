'''
Created on 2013-05-31

@author: jyeung
'''


import sys
import os
from scipy import stats
from utilities import read_write_gene_info, set_directories


# Folder names
inputdir = 'input'
outputdir = 'output'
tabledir = 'tables'

# Set directories
mydirs = set_directories.mydirs(inputdir, outputdir)

def standardize_rows(readwrite, rownames=True):
    NAcount = 0
    while True:
        try:
            row = readwrite.readnext()
        except StopIteration:
            print('%s rows read, %s rows written. Done.' %(readwrite.readrowcount, 
                                                           readwrite.writerowcount))
            print('%s rows contained NA, which were skipped.' %NAcount)
            break     
        if rownames == True:
            rowname = row[0]    # First element is rowname. 
            # if rowname is NA, skip it. 
            if rowname != 'NA':
                try:
                    gene_exprs = [float(i) for i in row[1:]]
                    gene_exprs_std = list(stats.zscore(gene_exprs))
                    readwrite.writenext([rowname] + gene_exprs_std)
                except ValueError:
                    # Assumes ValueError due to an NA value. 
                    NAcount += 1
                    print('ValueError, probably encountered NA value. Skipping.')
            else:
                print('Rowname is NA, skipping it.')
                NAcount += 1
                pass
        elif rownames == False:
            gene_exprs = [float(i) for i in row]
            gene_exprs_std = list(stats.zscore(gene_exprs))
            readwrite.writenext(gene_exprs_std)
    return None

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Input and output filename must be provided in the command line.')
        sys.exit()
    read_fname = sys.argv[1]
    write_fname = sys.argv[2]
    
    read_path = os.path.join(mydirs.outputdir, tabledir, read_fname)
    write_path = os.path.join(mydirs.outputdir, tabledir, write_fname)
    
    readwrite = read_write_gene_info.read_write_gene_info(read_path, write_path, header=False)
    
    with readwrite:
        # First two rows are headers, read them first before we standardize. 
        for _ in range(0, 2):
            header = readwrite.readnext()
            readwrite.writenext(header)
        # Now for each row, we standardize. 
        standardize_rows(readwrite, rownames=True)
    