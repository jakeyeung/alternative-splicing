'''
Created on 2014-01-18

@author: jyeung

Compares meme motifs with null gerp scores.
'''

import sys
from optparse import OptionParser
import pickle
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
from utilities import plot_functions, gerp_utilities

def get_dic_from_pklpath(pklpath):
    '''
    opens pickle file, loads into dic.
    '''
    with open(pklpath, 'rb') as pklfile:
        dic = pickle.load(pklfile)
    return dic

def get_gerp_scores(summary_dic, gerpkey='avg_rs_score'):
    '''
    From summary dic (output from summarize_meme_results_with_gerp_scores.py),
    retrieve all the gerp scores.
    
    Expect gerp score to be two subdics deep, key name "avg_rs_score"
    by default.
    
    I also expect the rs scores to be in list form (most of it is length 1, but
    it's not guaranteed)
    '''
    rs_scores = []
    for event in summary_dic:
        for region in summary_dic[event]:
            rs_scores += summary_dic[event][region][gerpkey]
    return rs_scores

def plot_distributions(meme_gerp_scores, null_gerp_scores, mylabels, mytitle):
    
    for gerp_scores, mylabel in zip([meme_gerp_scores, 
                                     null_gerp_scores], 
                                    mylabels):
        plot_functions.plot_density(gerp_scores, mytitle, mylabel)
    plt.legend()
    plt.show()
    
def main():
    usage = 'usage: %prog meme_gerp_genename_filepath output_filepath\n'\
        'Requires two input arguments:\n'\
        '1) pkl file from summarize_meme_results: non-null\n'\
        '2) pkl file from summarize_meme_results: null-mode\n'
    parser = OptionParser(usage=usage)   
    parser.add_option('-t', '--threshold', dest='score_threshold',
                      default=2.0,
                      help='Float, threshold for what one considers conserved.') 
    parser.add_option('-y', '--ymax', dest='ymax',
                      type='float',
                      default=0.03,
                      help='Y max for density plot')
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print 'Two arguments need to be specified in command line.\n'
        print usage
        sys.exit()
    non_null_pklpath = args[0]
    null_pklpath = args[1]
    # parse ops
    score_threshold = float(options.score_threshold)
    
    # get dics from pkl 
    non_null_dic = get_dic_from_pklpath(non_null_pklpath)
    null_dic = get_dic_from_pklpath(null_pklpath)
    
    non_null_gerp_scores = get_gerp_scores(non_null_dic, gerpkey='avg_rs_score')
    null_gerp_scores = get_gerp_scores(null_dic, gerpkey='avg_rs_score')
    
    plot_functions.plot_density([non_null_gerp_scores, null_gerp_scores], 
                                mytitle='Density plot of conservation scores', 
                                labels_lists=['MEME motifs', 'Controls'],
                                xlabel='GERP conservation score',
                                ylabel='Density',
                                xmin=-4, xmax=4,
                                ymax=options.ymax,
                                smoothness=0.15,
                                drawvline=score_threshold)
    
    # find how many conserved regions are in each.
    n_conserved_in_meme = \
        gerp_utilities.conserved_regions(non_null_gerp_scores, fraction=False, threshold=score_threshold)
    n_conserved_in_null = \
        gerp_utilities.conserved_regions(null_gerp_scores, fraction=False, threshold=score_threshold)
    n_total_in_meme = len(non_null_gerp_scores)
    n_total_in_null = len(null_gerp_scores)
    n_not_conserved_in_meme = n_total_in_meme - n_conserved_in_meme
    n_not_conserved_in_null = n_total_in_null - n_conserved_in_null
    
    print 'Threshold: %s' %score_threshold
    print 'Number of conserved elements: %s' %n_conserved_in_meme
    print 'Number of conserved elements found in control: %s' %n_conserved_in_null
    
    # Perform fisher's exact test
    oddsratio, pvalue = fisher_exact([[n_conserved_in_meme, 
                                       n_conserved_in_null], 
                                      [n_not_conserved_in_meme, 
                                       n_not_conserved_in_null]])
    print 'Fishers Exact Test, Oddsratio: %s. Pvalue: %s' %(oddsratio, pvalue)
    
    # plot distributions
    mylabels = ['Meme motifs', 'Control region']
    mytitle = 'Fraction of elements conserved compared to control region'
    # Plot bargraphs
    frac_conserved_meme = float(n_conserved_in_meme) / n_total_in_meme
    frac_conserved_null = float(n_conserved_in_null) / n_total_in_null
    myvals = [frac_conserved_meme, frac_conserved_null]
    plot_functions.plot_barplot(myvals, mytitle, mylabels, 
                                ylabel='Fraction of elements conserved', 
                                mytext1="%i/%i" \
                                    %(n_conserved_in_meme, 
                                      n_total_in_meme),
                                mytext2='%i/%i' %(n_conserved_in_null, 
                                                  n_total_in_null),
                                mytext3="*Fisher's Exact Test\nP-value=%.2e" %pvalue,
                                ymin=0,
                                ymax=1,
                                width=0.5)
    plt.show()
    
if __name__ == '__main__':
    main()