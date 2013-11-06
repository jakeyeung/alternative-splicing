'''
Created on 2013-09-18

@author: jyeung

Using existing functions to get gene names from gff3 file. 
'''

import sys
from append_gene_names_to_textfile import index_miso_annots, write_obj

def write_annot_dic_to_file(annot_dic, myobj):
    '''
    Given write_obj (csv.writer'd), take annotation dic and
    write each event into file as one column and gene name in another column.
    '''
    writecount = 0
    for event, genesymbol in annot_dic.iteritems():
        myobj.writeobj.writerow([event, genesymbol])
        writecount += 1
    print('%s rows written to file: %s' %(writecount, myobj.path))

def main(annot_file, output_file):
    annot_dic = index_miso_annots(annot_file)
    myobj = write_obj(output_file)
    write_annot_dic_to_file(annot_dic, myobj)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('MISO gff3 (version 2 and above) and output '\
              'file must be specified in command line.')
        sys.exit()
    annot_file = sys.argv[1]
    output_file = sys.argv[2]
    main(annot_file, output_file)