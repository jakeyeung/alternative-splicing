'''
Created on 2013-06-02

@author: jyeung
'''


import os


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
        
    def extract_first_col_all_files(self):
        '''
        Extract first column for each file in filelist.
        '''
        pass
        