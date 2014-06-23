'''
Created on 2014-02-22

@author: jyeung
'''

import sys
import csv
from optparse import OptionParser
from scipy.stats.stats import fisher_exact
from utilities.plot_functions import plot_barplot
import matplotlib.pyplot as plt

def count_anchor_results(results_file, colname='binding_regions'):
    '''
    Read results file, count rows containing
    binding_regions.

    Return rows containing binding regions and
    total number of rows counted.
    '''
    value_count = 0
    with open(results_file, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        header = myreader.next()
        for row_count, row in enumerate(myreader):
            value = row[header.index(colname)]
            if value is not '':
                value_count += 1
    return value_count, row_count

def main():
    usage = 'usage: %prog anchor_results.txt anchor_results_null.txt\n'\
        'Requires two input arguments:\n'\
        '1) Interesting anchor results, output from run_anchor_batch.py\n'\
        '2) Null anchor results, output from run_anchor_batch.py\n'
    parser = OptionParser(usage=usage)
    parser.add_option('-1', '--exon_label1', dest='exon_label1',
                      default='Exon label 1',
                      help='Exon label of anchor_results.txt.')
    parser.add_option('-2', '--exon_label2', dest='exon_label2',
                      default='Exon label 2',
                      help='Exon label of anchor_results_null.txt')
    parser.add_option('-t', '--title', dest='title',
                      default='Fraction of exons with predicted binding regions',
                      help='Title of plot.')
    (options, args) = parser.parse_args()
    if len(args) != 2:
        print 'Two arguments need to be specified in command line.\n'
        print usage
        sys.exit()
    anchor_results_path = args[0]
    anchor_results_null_path = args[1]
    exon_label1 = options.exon_label1
    exon_label2 = options.exon_label2
    mytitle = options.title

    # init dic with keys and empty lists
    anchor_dic = {}
    for key in ['binding', 'non_binding', 'total']:
        anchor_dic[key] = []

    for results in [anchor_results_path, anchor_results_null_path]:
        binding_count, total_count = count_anchor_results(results)
        non_binding_count = total_count - binding_count
        for key, val in zip(['binding', 'non_binding', 'total'],
                            [binding_count, non_binding_count, total_count]):
            anchor_dic[key].append(val)

    oddsratio, pvalue = \
        fisher_exact([anchor_dic['binding'], anchor_dic['non_binding']])

    print 'oddsratio: %s\npvalue: %s' %(oddsratio, pvalue)

    # plot distributions (from plot_meme_motif_null_comparison.py)
    mylabels = [exon_label1, exon_label2]
    # Plot bargraphs
    frac_binding = float(anchor_dic['binding'][0]) / anchor_dic['total'][0]
    frac_binding_null = float(anchor_dic['binding'][1]) / anchor_dic['total'][1]
    myvals = [frac_binding, frac_binding_null]
    plot_barplot(myvals, mytitle, mylabels,
                 ylabel='Fraction predicted binding regions',
                 mytext1="%i/%i" \
                    %(anchor_dic['binding'][0],
                      anchor_dic['total'][0]),
                  mytext2='%i/%i' %(anchor_dic['binding'][1],
                                    anchor_dic['total'][1]),
                  mytext3="*Fisher's Exact Test\nP-value=%.2e" %pvalue,
                  ymin=0,
                  ymax=1,
                  width=0.5)
    plt.show()

if __name__ == '__main__':
    main()
