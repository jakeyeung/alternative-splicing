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

def read_motifs_from_file(myfile, skipheader=True, rownames=False):
    '''
    Motif file should be tab delimited...
    
    Given a motif file (format is A C G U followed by
    a position-weighted matrix (PWM)),
    extract the fractions as a list, then return it.
    Skipheader: skips first line.
    Rownames = True -> removes first element in each row, since it is assumed
        to be a rowname. False does not remove any element in rows.
    '''
    mymotifs = []
    with open(myfile, 'rb') as readfile:
        # Skip header, contains A C G U, we know that alrdy.
        if skipheader == True:
            readfile.readline()
        if rownames == False:
            for l in readfile:
                mymotifs.append(l)
        elif rownames == True:
            for l in readfile:
                # Split by tab
                l_split = l.split('\t')
                # Re-join, but leaving out first element.
                l_norowname = '\t'.join(l_split[1:])
                mymotifs.append(l_norowname)
    return mymotifs