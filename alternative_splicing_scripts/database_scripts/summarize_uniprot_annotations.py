'''
Created on 2014-01-10

@author: jyeung

Parse annotation summary from annotate_genes_with_swissprot.py output
and plot barplots for each feature.

Also summarize total number of genes.

Inclusion / Exclusion is important, so we need to also read from
miso summary file.
'''

import sys
import csv
from optparse import OptionParser

def get_incl_excl_miso_events(miso_fpath, event_colname,
                              samp1_psi_colname, 
                              samp2_psi_colname):
    '''
    Read miso_fpath
    Requires input of samp1_psi_colname and samp2_psi_colname because
    the miso summary could be from bayes or T-test.
    event_colname also important since it differs between bayes and t-test
    
    Output a dictionary of form:
    {event_name: inclusion/exclusion}
    Inclusion/exclusion is with respect to sample 2.
    '''
    incl_str, excl_str = get_features_dic_keys()
    outdic = {}
    with open(miso_fpath, 'rb') as readfile:
        jreader = csv.reader(readfile, delimiter='\t')
        jheader = jreader.next()
        for row in jreader:
            event = row[jheader.index(event_colname)]
            samp1_psi = float(row[jheader.index(samp1_psi_colname)])
            samp2_psi = float(row[jheader.index(samp2_psi_colname)])
            # check if inclusion or exclusion by comparing samp2 to samp1.
            if samp2_psi > samp1_psi:
                # inclusion with respect to samp2
                if event not in outdic:
                    outdic[event] = incl_str
                else:
                    print 'Unexpected duplication in event: %s' %event
                    print 'Press enter to continue...'
                    raw_input()
            elif samp2_psi < samp1_psi:
                # exclusion
                if event not in outdic:
                    outdic[event] = excl_str
                else:
                    print 'Unexpected duplication in event: %s' %event
                    print 'Press enter to continue...'
                    raw_input()
            else:
                print 'Sample 2 PSI = Sample 1 PSI, '\
                    'press enter to label it inclusion.'
                print 'Samp2 PSI: %s\nSamp1 PSI: %s' %(samp2_psi, samp1_psi)
                raw_input()
    return outdic

def get_features_dic_keys():
    '''
    Return keys relevant for
    domains_dic
    '''
    incl_str = 'inclusion'
    excl_str = 'exclusion'
    return incl_str, excl_str

def count_n_feature(features_subdic, feature):
    '''
    Given feature and subdic, update feature 
    counts.
    Expect features_subdic to have form:
    {MOD_RES: 2}
    Increase count by 1 if matches feature.
    If subdic does not have feature as a subkey,
    initialize it with count at 1.
    '''
    if feature in features_subdic:
        features_subdic[feature] += 1
    else:
        features_subdic[feature] = 1
    return features_subdic

def get_features_dic(swissprot_annot_fpath, inclexcl_dic, 
                event_colname, feature_colname):
    '''
    Open swissprot annotation path, check each event
    if it is inclusion or exclusion, then put domain
    information in its respective dictionaries.
    Output dic format:
    {inclusion: [dom1, dom2, dom3...], exclusion: [dom4, dom5, dom6...]}
    '''
    features_dic = {}
    # Define keys in features dic
    incl_str, excl_str = get_features_dic_keys()
    # Init empty dics in features dic
    for key in [incl_str, excl_str]:
        features_dic[key] = {}
    
    with open(swissprot_annot_fpath, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        myheader = myreader.next()
        for row in myreader:
            event = row[myheader.index(event_colname)]
            feature = row[myheader.index(feature_colname)]
            # Check if inclusion or exclusion by looking up inclexcl_dic
            try:
                incl_or_excl = inclexcl_dic[event]
            except KeyError:
                print 'Could not find %s in incl_excl_dic' %event
                sys.exit()
            print incl_or_excl
            # put feature into features dic
            if incl_or_excl == incl_str:
                subdic = count_n_feature(features_dic[incl_str], feature)
                features_dic[incl_str].update(subdic)
            elif incl_or_excl == excl_str:
                subdic = count_n_feature(features_dic[excl_str], feature)
                features_dic[excl_str].update(subdic)
            else:
                print 'Expected %s to be %s or %s.' %(incl_or_excl, 
                                                      incl_str, excl_str)
    return features_dic

def main():
    usage = 'usage: %prog [opt] swissprot_annotated_results '\
        'miso_summary output_file'\
        '\nThree arguments must be specified in command line:\n'\
        '1) swissprot_annotated_results '\
        '(output from annotate_genes_with_swissprot.py)\n'\
        '2) miso_summary_file (bayes or T-Test).\n'\
        '3) Output file'
    parser = OptionParser(usage=usage)
    parser.add_option('--miso_event_colname', dest='miso_event_colname',
                      default='event_name',
                      help='Column name in miso summary file containing \n'\
                        'miso event. Could differ depending on whether '\
                        'file was from bayes or t-test.\n '\
                        'Default: "event_name". \nTry "event" '\
                        'if file was from t-test.')
    parser.add_option('--samp1_psi_colname', dest='samp1_psi_colname',
                      default='sample1_posterior_mean',
                      help='Sample 1 psi column name.\n'\
                        'Default: "sample1_posterior_mean"\n'\
                        'Could also be psi_mean1 if from t-test.')
    parser.add_option('--samp2_psi_colname', dest='samp2_psi_colname',
                      default='sample2_posterior_mean',
                      help='Sample 2 psi colname.\n'\
                        'Default: "sample2_posterior_mean"\n'\
                        'Could also be psi_mean2 if from t-test.')
    parser.add_option('--summary_event_colname', dest='summary_event_colname',
                      default='miso_event',
                      help='Swissprot annotated summary miso event colname.\n'\
                      'Default: "miso_event"\n')
    parser.add_option('--feature_colname', dest='feature_colname',
                      default='feature',
                      help='Swissprot annotated summary feature colname.\n'\
                      'Default: "feature"')
    (options, args) = parser.parse_args()
    if len(args) != 3:
        print 'Incorrect number of args specified.\n%s' %usage
        sys.exit()
    # Define positional args
    swissprot_annot_fpath = args[0]
    miso_summary_fpath = args[1]
    output_fpath = args[2]
    # Parse options
    miso_event_colname = options.miso_event_colname
    samp1_colname = options.samp1_psi_colname
    samp2_colname = options.samp2_psi_colname
    summary_event_colname = options.summary_event_colname
    feature_colname = options.feature_colname
    
    # Index miso summary, {misoevent: inclusion/exclusion}
    # allows separation of swissprot_annot_fpath into incl or excl.
    inclexcl_dic = get_incl_excl_miso_events(miso_summary_fpath, 
                                             miso_event_colname, 
                                             samp1_colname, samp2_colname)
    
    # domains_dic {inclusion: {domain: }, exclusion: {domain: }}
    features_dic = get_features_dic(swissprot_annot_fpath, inclexcl_dic, 
                                    summary_event_colname, feature_colname)
    print features_dic
    
if __name__ == '__main__':
    main()