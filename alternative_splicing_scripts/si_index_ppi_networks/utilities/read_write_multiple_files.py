'''
Created on 2013-06-12

@author: jyeung
'''


import sys
import csv


class read_write(object):
    '''
    Class for handling operations on multiple read files and one write file.
    '''


    def __init__(self, read_fullpath_list, write_fullpath):
        '''
        read_fullpath_list: list containing files to be read.
        write_fullpath: filepath to write.
        '''
        self.readpathlist = read_fullpath_list
        self.writepath = write_fullpath
        self.readrowcount = 0
        self.writerowcount = 0
        if self.writepath in self.readpathlist:
            print('Read and write path the same, exiting for your protection...')
            sys.exit()
            
    def __enter__(self):
        '''
        Open read and write paths as read and write objects.
        Initialize the reader objects and get reader colnames. 
        '''
        self.readfile_list = []
        for readpath in self.readpathlist:
            self.readfile_list.append(open(readpath, 'rb'))
        self.writefile = open(self.writepath, 'wb')
        
        self.reader_list = []
        for readfile in self.readfile_list:
            self.reader_list.append(csv.reader(readfile, delimiter='\t'))

        self.writer = csv.writer(self.writefile, delimiter='\t')
    
    def __exit__(self, exittype, exitvalue, exittraceback):
        '''
        Close the file.
        '''
        for readfile in self.readfile_list:
            readfile.close()
        self.writefile.close()
    
    def readnext(self):
        '''
        Read next row for ALL read files simultaneously. 
        '''
        row = []
        for reader in self.reader_list:
            row.append(reader.next())
        self.readrowcount += 1
        return row
    
    def writenext(self, row):
        self.writer.writerow(row)
        self.writerowcount += 1
        return row
        