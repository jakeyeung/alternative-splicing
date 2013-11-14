'''
Created on 2013-11-08

@author: jyeung
'''

import csv

def write_list_of_lists_to_file(llists, outputpath, header=False):
    '''
    Inputs: list of lists and output path.
    Consider each list within list of lists as a row, then write each
    row onto the output file path.
    
    If header != False, it is expected to be a list used to write header row.
    '''
    rowcount = 0
    with open(outputpath, 'wb') as writefile:
        jwriter = csv.writer(writefile, delimiter='\t')
        # Write header if it is not false.
        if header != False:
            jwriter.writerow(header)
        for l in llists:
            jwriter.writerow(l)
            rowcount += 1
    print '%s rows written to file: %s' %(rowcount, outputpath)
    return rowcount