'''
Created on 2014-06-16

@author: jyeung

Count number of unique genes in gff3

usage: python count_genes_in_gff3 mygff3file.gff3
'''

import sys
import csv


def extract_gsymbol(id_tag):
    '''
    Example input:
    ID=id;Name=name;gid=gid;ensg_id=eid;refseq_id=NA;gsymbol=IAH1
    
    Output:
    IAH1
    '''
    gsymbol = id_tag.split(';')[5]
    # split by equal sign
    gsymbol = gsymbol.split('=')[1]
    return gsymbol


def extract_miso_id(id_tag):
    '''
    Example input:
    ID=id;Name=name;gid=gid;ensg_id=eid;refseq_id=NA;gsymbol=IAH1
    
    Output:
    IAH1
    '''
    miso_id = id_tag.split(';')[0]
    # split by equal sign
    miso_id = id_tag.split('=')[1]
    return miso_id

    
def main():
    infile = sys.argv[1]
    
    gsymbol_dic = {}
    id_dic = {}
    
    n_genes = 0
    n_ids = 0
    
    with open(infile, 'rb') as readfile:
        jreader = csv.reader(readfile, delimiter='\t')
        jreader.next()  # skip first row
        for row in jreader:
            gene_or_exon = row[2]
            if gene_or_exon != 'gene':
                continue
            id_tag = row[8]
            # parse id tag
            gsymbol = extract_gsymbol(id_tag)
            miso_id = extract_miso_id(id_tag)
            if gsymbol not in gsymbol_dic:
                gsymbol_dic[gsymbol] = None
                n_genes += 1
            if id_tag not in id_dic:
                id_dic[miso_id] = None
                n_ids += 1
    
    print 'Number of genes in gff3 file: %s' % n_genes
    print 'Number of events in gff3 file: %s' % n_ids
            

if __name__ == '__main__':
    main()