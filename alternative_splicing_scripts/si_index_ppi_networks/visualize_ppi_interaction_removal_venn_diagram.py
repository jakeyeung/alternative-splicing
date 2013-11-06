'''
Created on 2013-07-08

@author: jyeung

After running remove_ppi_edges.py, we have a text file of 3 columns;
gene1, gene1, exist_or_not (1 or 0). 

We will plot a venn diagram between two tissue-specific ppi networks
 (e.g. adenocarcinoma vs neuroendocrine) and find the common and different
 protein interactions.
'''


import sys
import os
import csv
from utilities import set_directories, venn_diagram_plot


# Set constants
inputdir = 'input'
outputdir = 'output'

# Set directories
mydirs = set_directories.mydirs(inputdir, outputdir)

def find_diffs_two_ppi_networks(tissue_specific_ppi_fullpath1, 
                                tissue_specific_ppi_fullpath2, 
                                header=False):
    '''
    Read two tissue-specific PPI networks, which should contain three columns;
    gene1, gene2, and a binary value (0 or 1) indicating if the interaction exists 
    or not. Interaction/non-interaction derived by overlaying gene expression
    on ppi network.
    
    header by default is False.
    '''
    with open(tissue_specific_ppi_fullpath1, 'rb') as readfile1, \
    open(tissue_specific_ppi_fullpath2, 'rb') as readfile2:
        reader1 = csv.reader(readfile1, delimiter='\t')
        reader2 = csv.reader(readfile2, delimiter='\t')
        if header==True:
            reader1.next()
            reader2.next()
        elif header==False:
            pass
        else:
            sys.exit('header must be True or False: %s was used instead.' \
                     %header)
        
        # Initialize count to see if edges are same between two ppis
        same_count = 0
        different_count = 0
        rows_read = 0
        while True:
            try:
                row1 = reader1.next()
                row2 = reader2.next()
            except StopIteration:
                print('%s rows read, done reading ppi networks' %rows_read)
                break
            '''
            First column is gene1,
            Second column is gene2,
            Third column is binary (1 or 0) indicating if edge exists or not.
            We will check the two genes are identical in both row1 and row2,
            then see if edges are same.
            If edges are different (e.g. one is 1, other is 0), record as different.
            If edges are same (e.g. both are 1 or both are 0), record as same.
            '''
            if row1[0] == row2[0] and row1[1] == row2[1]:
                rows_read += 1
                if row1[2] == row2[2]:
                    same_count += 1
                else:
                    different_count += 1
            else:
                # Cannot compare if two genes are not the same.
                if row1[0] != row2[0]:
                    sys.exit('Each gene for each row in both ppi networks must be same.'\
                             '%s and %s are not the same.' %(row1[0], row2[0]))
                else:
                    sys.exit('Each gene for each row in both ppi networks must be same.'\
                             '%s and %s are not the same.' %(row1[1], row2[1]))
        return same_count, different_count


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Two tissue specific ppi networks filename (relative to output folder)'\
              'and tissue names must be specified in commandline.')
        sys.exit()
    tissue_specific_ppi_filename1 = sys.argv[1]
    tissue_specific_ppi_filename2 = sys.argv[2]
    tissue_name1 = sys.argv[3]
    tissue_name2 = sys.argv[4]
    
    # Construct full paths from relative filenames.
    # Assumes filenames are relative to output folder.
    tissue_specific_ppi_fullpath1 = os.path.join(mydirs.outputdir, 
                                                tissue_specific_ppi_filename1)
    tissue_specific_ppi_fullpath2 = os.path.join(mydirs.outputdir,
                                                 tissue_specific_ppi_filename2)
    
    same_count, diff_count = find_diffs_two_ppi_networks(tissue_specific_ppi_fullpath1, 
                                                               tissue_specific_ppi_fullpath2)
    venn_diagram_plot.plot_two_set_venn_diagram(diff_count, diff_count, 
                                                same_count, 
                                                [tissue_name1, tissue_name2])
    
    
    