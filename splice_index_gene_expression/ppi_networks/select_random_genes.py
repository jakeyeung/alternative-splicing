'''
Created on 2013-09-09

@author: jyeung

Given a list of genes, select a list of random ones.
Then we can run find_neighborhood_gene_exprs.py to find
differences.
'''


import sys
import csv
import random


def get_full_gene_list(gene_list_filename):
    # Define constants and lists
    gene_colname = 'gene'
    gene_list = []
    
    with open(gene_list_filename, 'rb') as myfile:
        myreader = csv.reader(myfile, delimiter='\t')
        header = myreader.next()
        for row in myreader:
            gene = row[header.index(gene_colname)]
            gene_list.append(gene)
    print('Loaded %s genes into list from file: %s' %(len(gene_list), 
                                                   gene_list_filename))
    return gene_list

def choose_from_list(mylist, n):
    random.seed(1)    # Set seed
    return random.sample(mylist, n)

def write_list_to_file(mylist, outputpath):
    with open(outputpath, 'wb') as writefile:
        writecount = 0
        for i in mylist:
            writefile.write('%s\n' %i)
            writecount += 1
    print('Written %s lines to file in %s' %(writecount, outputpath))
    
def main(gene_list_filename, output_filename):
    gene_list = get_full_gene_list(gene_list_filename)
    random_genes = choose_from_list(gene_list, n=450)
    write_list_to_file(random_genes, output_filename)
        
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('List of genes and output filepath '\
              'must be specified in command line.')
        sys.exit()
    gene_list_filename = sys.argv[1]
    output_filename = sys.argv[2]
    main(gene_list_filename, output_filename)