'''
Created on 2013-10-23

@author: jyeung
Get background list ready for IPA
'''

import sys 
import csv

def row_is_gene(row):
    '''
    Check if row is gene in third column.
    '''
    gene_index = 2
    if row[gene_index] == 'gene':
        return True
    else:
        return False

def get_gff3_genes(gff3_file):
    '''
    Get gene list.
    '''
    # Def constants
    info_index = 8
    
    gene_list = []
    with open(gff3_file, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        myreader.next()    # Skip header
        for row in myreader:
            if row_is_gene(row):
                info = row[info_index]
                # Split by semicolon and equal sign.
                colsplit = info.split(';')[5]
                genename = colsplit.split('=')[1]
                if genename != 'NA':
                    if genename not in gene_list:
                        gene_list.append(genename)
            else:
                pass
    return gene_list

def get_column_from_file(myfile, col_index, header=True):
    '''
    Read file, col_index.
    '''
    mylist = []
    with open(myfile, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        if header==True:
            myreader.next()    # Skip header
        for row in myreader:
            mylist.append(row[col_index])
    return mylist

def mark_genes_in_list(bg_list, ds_list, output_file):
    '''
    If DS gene in bg_list, mark it as 'X'.
    Write to output file.
    '''
    mark_count = 0
    with open(output_file, 'wb') as writefile:
        mywriter = csv.writer(writefile, delimiter='\t')
        for bg_gene in bg_list:
            if bg_gene in ds_list:
                mywriter.writerow([bg_gene, 'X'])
                mark_count += 1
            else:
                mywriter.writerow([bg_gene, ''])
    return mark_count

def main():
    gff3_file = sys.argv[1]
    ds_genes_file = sys.argv[2]
    output_file = sys.argv[3]
    
    # Get bg genes
    bg_list = get_gff3_genes(gff3_file)
    print '%s genes indexed.' %len(bg_list)
    
    # Get DS genes
    ds_list = get_column_from_file(ds_genes_file, 1)
    
    mark_count = mark_genes_in_list(bg_list, ds_list, output_file)
    print '%s genes marked as DS in BG list.' %mark_count

if __name__ == '__main__':
    main()