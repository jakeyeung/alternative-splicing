'''
Created on 2013-05-30

@author: jyeung
From a large table containing splice index, probe expression gene expression,
collect ONLY gene expression, set it up into an output text table to be used 
for OptDis. 
'''


import os
import sys
from utilities import set_directories, read_write_gene_info


# Set constants

# Folder names
inputdir = 'input'
cohortdir = 'takeda'
outputdir = 'output'
tabledir = 'tables'

# Input column names
mean_gene_exprs_colname = 'mean'
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

# full_sample_list containing all samples, for reference only. 
'''
full_sample_list = ['X946_Urethra', 
                    'X972.2.Penila',
                    'AB352',
                    'C42.RNA',
                    'LN.AI.Luc',
                    'X963.Lpa.JP',
                    'X963.L.LN',
                    'X1005',
                    'X890L',
                    'X890.LN',
                    'X945',
                    'X945L.LN',
                    'X961']
'''

# Set directories
mydirs = set_directories.mydirs(inputdir, outputdir)


def extract_and_write_gene_exprs(readwrite):    
    # Read rows, extract gene information. 
    with readwrite:
        # Get index of gene exprs start and gene exprs end.
        # The end index, when used with row[start:end] syntax
        # gives always the column name BEFORE mean_gene_exprs_colname, which 
        # is exactly what we want. 
        exprs_start_index = readwrite.inputcolnames.index(mean_gene_exprs_colname)-full_sample_length
        exprs_end_index = readwrite.inputcolnames.index(mean_gene_exprs_colname)
        rownames_subset = readwrite.inputcolnames[exprs_start_index:exprs_end_index]
        
        # Write headers
        # We write not only sample names, but we write the class the sample is in. 
        # Write sample names
        # writecols = [gene_label_colname]
        classnames = ['']
        samplenames = ['']
        for tumor_class, sample_name in classifier_dict.iteritems():
            for s in sample_name:
                classnames.append(tumor_class)
                samplenames.append(s)
        readwrite.writenext(samplenames)
        readwrite.writenext(classnames)
        
        '''
        Read and write gene expression information.
        '''
        count = 1
        writecount = 0
        prev_gene_name = 0
        while True:
            try:
                row = readwrite.readnext()
                gene_name = row[readwrite.inputcolnames.index(gene_symbol_colname)]
                # subset from mean_gene_exprs_colname minus full sample length
                #  to mean_gene_exprs_colname
                rowsubset = row[exprs_start_index:exprs_end_index]
            except StopIteration:
                print('%s rows read, %s rows written. Done.' %(readwrite.readrowcount, writecount))
                break
            if gene_name != prev_gene_name:
                # Extract only gene expressions in sample list.
                writecount += 1
                gene_exprs_list = [gene_name]
                for s in sample_list:
                    gene_exprs = rowsubset[rownames_subset.index(s)]
                    gene_exprs_list.append(gene_exprs)
                readwrite.writenext(gene_exprs_list)
            prev_gene_name = gene_name
            count += 1
    return None

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Input and output filename must be provided in the command line.')
        sys.exit()
    read_fname = sys.argv[1]
    write_fname = sys.argv[2]
    
    read_path = os.path.join(mydirs.inputdir, cohortdir, read_fname)
    write_path = os.path.join(mydirs.outputdir, tabledir, write_fname)
    
    readwrite = read_write_gene_info.read_write_gene_info(read_path, write_path)
    
    extract_and_write_gene_exprs(readwrite)
    

        
    
    