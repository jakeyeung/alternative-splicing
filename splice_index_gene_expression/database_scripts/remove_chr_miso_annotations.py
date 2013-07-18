'''
Created on 2013-07-16

@author: jyeung

Remove 'chr' from chromosome indicators in the miso database.
Why? It seems Mark Rubin's cohort was aligned in such a way so the chromosomes 
were labelled 1, 2, 3, 4, 5... 21, 22, X, Y, MT. 
'''


import sys
import csv


def remove_chr_prefix(input_file, output_file, first_col_only=True):
    '''
    Take input_file and finds string in first column and removes 'chr' from that 
    string.
    
    Assumes input_file has a header in first line.
    
    output_file prints same as input_file but first column has chr removed.
    
    If first_col_only == True, then it will only remove 'chr' from the first column.
    If False, it will also remove from the 9th column, which indicates an id. 
    '''
    with open(input_file, 'rb') as read_file, \
    open(output_file, 'wb') as write_file:
        # Write first line to write_file
        write_file.write(read_file.readline())
        # Read rest of file as csv reader
        reader = csv.reader(read_file, delimiter='\t')
        writer = csv.writer(write_file, delimiter='\t')
        rowcount = 0
        for row in reader:
            '''
            Assume first column contains chromosome values and it is
            prefixed by 'chr'. Replace 'chr' with '' then write
            the modified row to file. 
            
            Also replace 'chr' with '' for index 8 containing ID and Name information.
            '''
            row[0] = row[0].replace('chr', '')
            if first_col_only == False:
                row[8] = row[8].replace('chr', '')
            elif first_col_only == True:
                pass
            else:
                sys.exit('%s: first_col_only must be either False or True.')
            writer.writerow(row)
            rowcount += 1
        print('%s rows iterated.' %rowcount)


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Alternatively spliced database .gff file '\
              'and output write file must be provided in command line.')
        sys.exit()
    read_file = sys.argv[1]
    write_file = sys.argv[2]
    
    remove_chr_prefix(read_file, write_file)
    
    