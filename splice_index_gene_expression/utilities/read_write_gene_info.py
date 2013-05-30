'''
Created on 2013-05-30

@author: jyeung
'''


import csv


class read_write_gene_info(object):
    '''
    For reading and writing gene information.
    '''
    
    def __init__(self, read_fullpath, write_fullpath):
        '''
        Constructor
        '''
        self.readpath = read_fullpath
        self.writepath = write_fullpath
        self.readrowcount = 0
        self.writerowcount = 0
        
    def __enter__(self):
        '''
        Open read and write paths as read and write objects.
        Initialize the reader objects and get reader colnames. 
        '''
        self.readfile = open(self.readpath, 'rb')
        self.writefile = open(self.writepath, 'wb')
        
        self.reader = csv.reader(self.readfile, delimiter='\t')
        self.writer = csv.writer(self.writefile, delimiter='\t')
        
        self.inputcolnames = self.reader.next()
        # Initialize output header in writecolnames() function.
        
    def __exit__(self, exittype, exitvalue, exittraceback):
        '''
        Close the file.
        '''
        self.readfile.close()
        self.writefile.close()
        
    def writecolnames(self, *colnames):
        header = []
        for cname in colnames:
            header.append(cname)
        self.outputheaders = header
        self.writer.writerow(header)
    
    def readnext(self):
        row = self.reader.next()
        self.readrowcount += 1
        return row
    
    def writenext(self, row):
        self.writer.writerow(row)
        self.writerowcount += 1
        return row
    
        
    
    
        