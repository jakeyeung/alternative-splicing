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


def remove_gene_and_mrna_from_gff3(input_file, output_file):
    '''
    From input file, read each line, if it contains gene
    or mRNA in the third column, skip. 
    Otherwise (if it contains exon), write the file to 
    the output.
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
        for row in reader:
            rowcount += 1
            if row[2] == 'gene' or row[2] == 'MRNA':
                pass
            elif row[2] == 'exon':
                writer.writerow(row)
                writecount += 1
        print('%s rows read in %s. %s rows written to %s' %(rowcount, 
                                                            input_file, 
                                                            writecount, 
                                                            output_file))


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Alternatively spliced database .gff file '\
              'and output write file must be provided in command line.')
        sys.exit()
    read_file = sys.argv[1]
    write_file = sys.argv[2]
    
    remove_gene_and_mrna_from_gff3(read_file, write_file)
    
    