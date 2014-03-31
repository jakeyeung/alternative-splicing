'''
Created on 2014-01-15

@author: jyeung

Quick script: convert ensembl gene id to gene symbol.
'''

import sys
import csv

def main():
    exprs_file = sys.argv[1]
    ensembl_gene_symbol_file = sys.argv[2]
    out_file = sys.argv[3]
    
    # index ensembl gene symbol
    ensembl_gene_dic = {}
    with open(ensembl_gene_symbol_file, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        # first column: ensemblID
        # second column: gene symbol
        for row in myreader:
            ensemblid = row[0]
            genename = row[1]
            ensembl_gene_dic[ensemblid] = genename
    
    # open write obj
    writefile = open(out_file, 'wb')
    mywriter = csv.writer(writefile, delimiter='\t')
    
    # open file, r
    with open(exprs_file, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter=' ')
        header = myreader.next()
        # Remove quotes from headers
        header = [i.strip('"') for i in header]
        # write header to file
        mywriter.writerow(header)
        # first column: ensemblID
        writecount = 0
        for row in myreader:
            ensemblid = row[0].strip('"')
            try:
                genename = ensembl_gene_dic[ensemblid]
                row[0] = genename
                mywriter.writerow(row)
                writecount += 1
            except KeyError:
                pass
    writefile.close()
    print 'Written %s rows to file: %s' %(writecount, out_file)

if __name__ == '__main__':
    main()