'''
Created on 2014-03-14

@author: jyeung
Given motif summary (html) file,
grab all the event ids.

Read miso output file (with gene name appended)
to correlate event ids with gene names.

Write output
'''

import sys
import csv
from optparse import OptionParser
from check_motif_cooperativity import get_seq_names_from_meme_results

def index_gene_names(miso_path, 
                     event_colname='event_name', 
                     gene_colname='gsymbol'):
    '''
    Input: miso_path containing seq ids and gene names
    Output: 
        Dictionary of form:
        {seq_id: mygene}
    '''
    outdic = {}
    with open(miso_path, 'rb') as readfile:
        jreader = csv.reader(readfile, delimiter='\t')
        header = jreader.next()
        for row in jreader:
            # get event and gene names
            event = row[header.index(event_colname)]
            gene = row[header.index(gene_colname)] 
            # put them into dictionary
            outdic[event] = gene
    return outdic

def main():
    usage = 'usage: %prog meme_results_html motif_number '\
        'miso_gene_name_file outputfile'\
        '\n'\
        'Four args must be specified in commandline: \n'\
        '1) Meme HTML file containing meme results.\n'\
        '2) Integer indicating motif number to match.'\
        ' Must not exceed number of motifs in meme result.\n'\
        '3) MISO file with gene name appended, '\
            'linking miso event with genename\n'\
        '4) Output file\n'\
        '-h to display this help message.\n'
    parser = OptionParser(usage=usage)
    parser.add_option('-e', '--event_colname', dest='event_colname', 
                      default='event_name',
                      help='Event colname in miso file, default "event_name"')
    parser.add_option('-g', '--gene_colname', dest='gene_colname', 
                      default='gsymbol',
                      help='Gene symbol colname in miso file, '\
                        'default "gsymbol"')
    (options, args) = parser.parse_args()
    if len(args) != 4:
        print 'Incorrect number of parameters specified.'
        print usage
        sys.exit()
    meme_html_path = args[0]
    motif_number = args[1]
    miso_path = args[2]
    output_path = args[3]
    
    event_gene_dic = index_gene_names(miso_path, 
                                      options.event_colname, 
                                      options.gene_colname)
    print '%s events indexed into dic' %len(event_gene_dic.keys())
    
    seq_ids = get_seq_names_from_meme_results(meme_html_path, motif_number)
    print '%s sequences matched to motifs %s' %(len(seq_ids), motif_number)
    
    matched_genes = []
    for seq in seq_ids:
        matched_genes.append(event_gene_dic[seq])
    
    print 'Sequences and gene names outputted to: %s' %output_path
    
    with open(output_path, 'wb') as outfile:
        jwriter = csv.writer(outfile, delimiter='\t')
        # write header: seq, genename
        jwriter.writerow([options.event_colname, options.gene_colname])
        for count, seq in enumerate(seq_ids):
            jwriter.writerow([seq, event_gene_dic[seq]])
    print '%s rows written to: %s' %(count+1, output_path)

if __name__ == '__main__':
    main()