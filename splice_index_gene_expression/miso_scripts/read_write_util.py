'''
Created on 2013-07-26

@author: jyeung
'''

import csv
import sys


class read_write(object):
    '''
    For reading and writing gene information.
    '''
    
    def __init__(self, read_fullpath, write_fullpath, header=True):
        '''
        Constructor
        '''
        self.readpath = read_fullpath
        self.writepath = write_fullpath
        self.readrowcount = 0
        self.writerowcount = 0
        self.header = header
        if self.readpath == self.writepath:
            print('Read and write path the same, exiting for your protection...')
            sys.exit()
        
    def __enter__(self):
        '''
        Open read and write paths as read and write objects.
        Initialize the reader objects and get reader colnames. 
        '''
        self.readfile = open(self.readpath, 'rb')
        self.writefile = open(self.writepath, 'wb')
        
        self.reader = csv.reader(self.readfile, delimiter='\t')
        self.writer = csv.writer(self.writefile, delimiter='\t')
        
        if self.header == True:
            self.inputcolnames = self.reader.next()
        elif self.header == False:
            pass
        else:
            print('Error, header arg accepts only True or False.')
            sys.exit()
        # Initialize output header in writecolnames() function.
        
    def __exit__(self, exittype, exitvalue, exittraceback):
        '''
        Close the file.
        '''
        self.readfile.close()
        self.writefile.close()
    
    def readnext(self):
        row = self.reader.next()
        self.readrowcount += 1
        return row
    
    def writenext(self, row):
        self.writer.writerow(row)
        self.writerowcount += 1
        return row