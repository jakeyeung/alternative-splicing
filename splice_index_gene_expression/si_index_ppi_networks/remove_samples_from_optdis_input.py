'''
Created on 2013-06-11

@author: jyeung
Reads an optdis input and writes to file the same
optdis input but with certain samples removed. 
'''


import sys
import os
from utilities import set_directories, read_write_gene_info


# Set constants
inputdir = 'input'
outputdir = 'output'

# Set directories
mydirs = set_directories.mydirs(inputdir, outputdir)

def remove_sample_names_from_optdis_output(optdis_fullpath, 
                                           output_fullpath, 
                                           samples_to_remove):
    '''
    Read optdis inputs, remove user-defined samples, output to file. 
    Inputs:
        optdis_fullpath: optdis inputs, contains TWO headers, first header is
            sample names, second header is the class that sample is in. Length
            of each row should be the same (same number of elements), firs two
            headers may contain a blank first element, but lengths should still
            be same as the rows in the body. 
        output_fullpath: path to write output file.
        samples_to_remove: list of sample names to be removed from optdis inputs. 
    '''
    optdis_readwrite = read_write_gene_info.read_write_gene_info(optdis_fullpath, 
                                                                 output_fullpath, 
                                                                 header=False)
    with optdis_readwrite:
        '''
        Find which column to delete based on samples_to_remove and first row. 
        '''
        print('Reading optdis inputs and writing to %s, removing %s for each row' \
              %(output_fullpath, samples_to_remove))
        samplename_list = optdis_readwrite.readnext()    # first row
        
        # Get indices to remove
        indices_to_remove = []
        for s in samples_to_remove:
            indices_to_remove.append(samplename_list.index(s))
        
        indices_to_remove = sorted(indices_to_remove, reverse=True)
        
        # Remove columns by index.
        for i in indices_to_remove:
            samplename_list.pop(i)
        
        optdis_readwrite.writenext(samplename_list)
        
        # Now do the same for all other rows. 
        while True:
            try:
                row = optdis_readwrite.readnext()
            except StopIteration:
                print('%s rows read, %s rows written, breaking...' \
                      %(optdis_readwrite.readrowcount, 
                        optdis_readwrite.writerowcount))
                break
            for i in indices_to_remove:
                row.pop(i)
            optdis_readwrite.writenext(row)
        
        return None

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('OptDis inputs, output filename, comma-separated samplenames must be specified in command line.'\
              '\nFilenames are relative to output.')
        sys.exit()
    
    optdis_partialpath = sys.argv[1]
    output_partialpath = sys.argv[2]
    sample_names = sys.argv[3]
    
    optdis_fullpath = os.path.join(mydirs.outputdir, optdis_partialpath)
    output_fullpath = os.path.join(mydirs.outputdir, output_partialpath)
    
    # Parse comma-separated samples.
    samples_to_remove = sample_names.split(',')
    
    # Read optdis file, remove unwanted samples, write to file.
    remove_sample_names_from_optdis_output(optdis_fullpath, 
                                           output_fullpath, 
                                           samples_to_remove)