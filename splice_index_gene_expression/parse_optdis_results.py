'''
Created on 2013-06-02

@author: jyeung
'''


import sys
import os
from utilities import set_directories, files_in_directory


# Set constants
inputdir = 'input'
outputdir = 'output'

# Set directories
mydirs = set_directories.mydirs(inputdir, outputdir)


def get_files_in_directory(module_fullpath):
    '''
    Get a list of files in a directory.
    '''
    filelist = os.listdir(module_fullpath)
    return filelist
    
    

def get_genes_from_filelist(filelist):
    '''
    The modules directory contains n number of files corresponding to
    n number of subnetworks. 
    
    Each subnetwork contains between 3 to 8?? number of genes. 
    
    The point here is to create a textfile that consolidates all this
    information so it will contain TWO COLUMNS:
        - a list of genes
        - to which module the gene belongs
    '''
    pass
    
    
    

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Module directory must be provided in command line.')
        sys.exit()
    module_directory = sys.argv[1]
    
    module_fullpath = os.path.join(mydirs.outputdir, module_directory)
    
    optdis_results = files_in_directory.files_in_directory(module_fullpath)
    
    print optdis_results.dir
    print optdis_results.filelist
    
    
    