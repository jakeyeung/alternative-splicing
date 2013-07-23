'''
Created on 2013-07-16

@author: jyeung

The gff3 downloaded from miso annotations
(http://genes.mit.edu/burgelab/miso/docs/annotation.html)
contains, for each alternative splice event, a gene,
two mRNA, and a list of their corresponding exons.

We want to remove the gene and two mRNAs from the gff3 so that
when we view the gff3 on IGV, it would be viewable as separate exons.
'''


import sys
import csv


def remove_gene_and_mrna_from_gff3(input_file, output_file, column_index, 
                                   list_to_remove, list_to_keep):
    '''
    Usage:
        Reads each row from input_file (gff3), at column_index of each row,
        make a decision whether to write the row to output depending on if it
        contains certain strings. 
    Inputs:
    input_file: gff3 file from which we want to remove some rows.
    output_file: gff3 file with certain rows removed.
    column_index: which column to check for a certain string. If it contains
        a string, remove row if that string is in list_to_remove. 
        keep row if that string is in list_to_keep.
    list_to_remove: a list of strings which will be searched at row[column_index].
        if it contains one of these strings, we will remove that row from the 
        output (remove meaning we will not write it to output file).
    list_to_keep: a list of strings to be searched at row[column_index].
        if it contains one of these strings, we will write that row to output.
    Output:
        an output file containing a subset of input_file rows. 
    '''
    with open(input_file, 'rb') as read_file, \
    open(output_file, 'wb') as write_file:
        # Write first line to write_file
        write_file.write(read_file.readline())
        # Read rest of file as csv reader
        reader = csv.reader(read_file, delimiter='\t')
        writer = csv.writer(write_file, delimiter='\t')
        rowcount = 0
        writecount = 0
        warningcount = 0
        for row in reader:
            rowcount += 1
            if row[column_index] in list_to_remove:
                pass
            elif row[column_index] in list_to_keep:
                writer.writerow(row)
                writecount += 1
            else:
                print('Warning: %s neither gene, mRNA or exon.' %row[2])
                print('Continuing anyway...')
                warningcount += 1
                pass
        print('%s rows read in %s. %s rows written to %s' %(rowcount, 
                                                            input_file, 
                                                            writecount, 
                                                            output_file))
        print('%s times a warning was triggered.' %warningcount)


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('1) Alternatively spliced database .gff file, '\
              '2) output write file, '\
              '3) column index in which to search, '\
              '4) comma-separated strings to search to filter out, '\
              '5) comma-separated strings to search to keep, must be provided in command line.')
        sys.exit()
    read_file = sys.argv[1]
    write_file = sys.argv[2]
    try:
        column_index = int(sys.argv[3])
    except ValueError:
        sys.exit('Unable to convert %s to an integer' %sys.argv[3])
    csv_strings_for_filter = sys.argv[4]
    csv_strings_for_keeping = sys.argv[5]
    
    list_to_remove = csv_strings_for_filter.split(',')
    list_to_keep = csv_strings_for_keeping.split(',')
    
    remove_gene_and_mrna_from_gff3(read_file, write_file, 
                                   column_index, list_to_remove, 
                                   list_to_keep)