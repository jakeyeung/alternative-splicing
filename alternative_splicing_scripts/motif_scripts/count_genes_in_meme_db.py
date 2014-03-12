'''
Created on 2014-01-19

@author: jyeung

Read custom meme database and return how many
unique genes are in the database.
'''

import sys

def main():
    if len(sys.argv) != 2:
        print 'Meme database filepath must be input in commandline.'
        sys.exit()
    meme_db = sys.argv[1]
        
    genes = []
    
    with open(meme_db, 'rb') as readfile:
        for line in readfile:
            if line.startswith('MOTIF'):
                # Example of line starting with motif:
                # MOTIF MBNL3,M320_0.6,I
                motif_id = line.split(' ')[1]    # 2nd element
                gene = motif_id.split(',')[0]
                genes.append(gene)
    genes = list(set(genes))
    
    print '%s genes in meme database.' %len(genes)
    print genes

if __name__ == '__main__':
    main()