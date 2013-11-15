'''
Created on 2013-11-14

@author: jyeung
'''

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

def read_motifs_from_file(myfile, skipheader=True):
    '''
    Given a motif file (format is A C G U followed by
    a position-weighted matrix (PWM)),
    extract the fractions as a list, then return it.
    '''
    mymotifs = []
    with open(myfile, 'rb') as readfile:
        # Skip header, contains A C G U, we know that alrdy.
        if skipheader == True:
            readfile.readline()
        for l in readfile:
            mymotifs.append(l)
    return mymotifs