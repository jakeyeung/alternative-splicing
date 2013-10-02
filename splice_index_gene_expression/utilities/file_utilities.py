'''
Created on 2013-09-27

@author: jyeung
'''

import csv
import sys

class csv_obj(object):
    '''
    Prepare read obj for reading or writing file
    '''
    def __init__(self, filename, filetype):
        self.filepath = filename
        self.type = filetype
        if self.type != 'read' and self.type != 'write':
            print('filetype must be either "read" or'\
                  ' "write". %s found...' %filetype)
            sys.exit()
        
    def __enter__(self):
        if self.type == 'read':
            self.file = open(self.filepath, 'rb')
            self.readobj = csv.reader(self.file, delimiter='\t')
        elif self.type == 'write':
            self.file = open(self.filepath, 'wb')
            self.writeobj = csv.writer(self.file, delimiter='\t')
        else:
            print('filetype must be either '\
                  '"read" or "write".. %s found.' %self.type)
            sys.exit()
    
    def __exit__(self, jtype, jvalue, tb):
        self.file.close()
        

def get_samples_from_file(sample_path):
    '''
    From a text file of one column containing sample names, 
    iterate rows and get the samples.
    '''
    samples = []
    with open(sample_path, 'rb') as readfile:
        reader = csv.reader(readfile, delimiter='\t')
        for row in reader:
            samples += row
    return samples