'''
Created on 2013-06-02

@author: jyeung
'''


import os
import csv


class files_in_directory(object):
    '''
    Class for handling multiple files in a directory.
    '''


    def __init__(self, directory):
        '''
        Constructor
        '''
        self.dir = directory
        self.filelist = os.listdir(directory)
        
    def extract_first_col_all_files(self, sortresults=False):
        '''
        Extract first column for each file in filelist.
        '''
        gene_list = []
        for f in self.filelist:
            with open(os.path.join(self.dir, f), 'rb') as file_obj:
                filereader = csv.reader(file_obj, delimiter='\t')
                for row in filereader:
                    genename = row[0]
                    module_number = f
                    gene_list.append((genename, module_number))
        if sortresults == False:
            self.gene_list = gene_list
        elif sortresults == True:
            self.gene_list = sorted(gene_list, key=lambda x: int(x[1]))
        return gene_list
    
    def write_genes_to_file(self, output_fullpath):
        with open(output_fullpath, 'wb') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            for gene in self.gene_list:
                genename = gene[0]
                modulenumber = gene[1]
                writer.writerow([genename, modulenumber])
        return None
                