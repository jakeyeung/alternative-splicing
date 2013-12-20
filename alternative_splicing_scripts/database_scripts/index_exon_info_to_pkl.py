'''
Created on 2013-12-19

@author: jyeung

Code adotped from cbayly
'''

import sys 
from convert_dna_fasta_to_protein_fasta import make_dictionary

def main():
    ensdatabase = sys.argv[1]
    ensdictionary = sys.argv[2]
    
    if len(sys.argv) < 3:
        print '3 arguments must be specified in command line.'
        sys.exit()
    
    print 'Making dictionary from %s' %ensdatabase
    _, count = make_dictionary(ensdatabase, ensdictionary)
    print 'Dictionary with %s entries saved to %s' %(count, ensdictionary)
    
if __name__ == '__main__':
    main()
