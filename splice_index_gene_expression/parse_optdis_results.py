'''
Created on 2013-06-02

@author: jyeung
'''


import sys
import os
import csv
from utilities import set_directories, files_in_directory


# Set constants
inputdir = 'input'
outputdir = 'output'
gene_column_index = 0    # First column of row in output contains gene names. 

# Set directories
mydirs = set_directories.mydirs(inputdir, outputdir)


def get_files_in_directory(module_fullpath):
    '''
    Get a list of files in a directory.
    '''
    filelist = os.listdir(module_fullpath)
    print filelist
    return filelist
    
def get_genes_from_file(filepath, initial_gene_list):
    '''
    From a filepath, read its file and return its genes.
    
    initial_gene_list contains a starter list with which new
    genes will be appended to the end. This makes this piece of code
    loopable. 
    ''' 
    with open(filepath, 'rb') as file_obj:
        filereader = csv.reader(file_obj, delimiter='\t')
        for row in filereader:
            initial_gene_list.append(row[0])
    return initial_gene_list

def get_genes_from_filelist(filelist):
    '''
    Input: filelist is a list of filenames to be read and information
    extracted. 
    
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
    output_dir = os.path.dirname(module_fullpath)
    output_filename = sys.argv[2]
    output_fullpath = os.path.join(output_dir, output_filename)
    optdis_results = files_in_directory.files_in_directory(module_fullpath)
    
    gene_list = optdis_results.extract_first_col_all_files(sortresults=True)
    
    optdis_results.write_genes_to_file(output_fullpath)
    
    print('%s genes written to file.' %len(gene_list))
    
    '''
    outputpath = os.path.join(module_fullpath, 'full_list.txt')
    
    with open(outputpath, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        for gene in gene_list:
            writer.writerow([gene])
    '''
    
    