'''
Created on 2013-11-14

@author: jyeung
'''

import sys
import csv

def store_textfile_as_list(textfile, collapse_list=True):
    '''
    Given a textfile, take each row as an element, append it to a list. 
    Collapse list = True: useful if textfile only contains one column per row.
    It then returns a list instead of list of lists.
    Return the list.
    
    '''
    mylist = []
    with open(textfile, 'rb') as myfile:
        myreader = csv.reader(myfile, delimiter='\t')
        for row in myreader:
            mylist.append(row)
    if collapse_list == True:
        mylist = [item for sublist in mylist for item in sublist]
    elif collapse_list == False:
        pass
    else:
        print 'Collapse_list must be either True or False (boolean)'
    return mylist