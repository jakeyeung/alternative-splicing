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
from utilities import plot_utils
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

def prepare_dic_for_barplot(incl_excl_dic, ignore_list=[]):
    '''
    From dictionary made from
    summarize_uniprot_annotations.py,
    prepare for 0plot barchart to compare inclusion
    vs exclusion.
    
    input dic is of form:
    {inclusion: {feature1: n1},
    exclusion: {feature2: n2}}.
    
    We want vectors of just features (for x-ticks)
    and values in each featuer (y-values)
    
    ignore_list: list of features to ignore.
    '''
    # Get full list of features, make it non-redundant
    features_list = []
    incl_str, excl_str = get_features_dic_keys()
    incl_excl_list = [incl_str, excl_str]
    
    for incl_or_excl in incl_excl_dic:
        features_list += incl_excl_dic[incl_or_excl].keys()
    # Remove redundancies.
    features_list = list(set(features_list))
    
    # Output features list may be smaller than features_list
    # due to ignore_list, so init the output features list
    out_features_list = []
    
    # Build bar heights
    n_features1 = []
    n_features2 = []
    for incl_excl, features_vector in zip(incl_excl_list, 
                                          [n_features1, 
                                           n_features2]):
        for feature in features_list:
            # ignore features in ignore_list
            if feature not in ignore_list:
                out_features_list.append(feature)
                try:
                    features_vector.append(incl_excl_dic[incl_excl][feature])
                except KeyError:
                    # key may not be in dic, then just make it 0
                    features_vector.append(0)
            else:
                print 'Ignoring feature: %s' %feature
    return out_features_list, n_features1, n_features2

def get_exon_count(protein_summary_file):
    '''
    Read protein summary file
    (e.g. exon_2_exclusion.protein.summary).
    This is an output from create_dna_protein_summary_file.py
    '''
    with open(protein_summary_file, 'rb') as readfile:
        jreader = csv.reader(readfile, delimiter='\t')
        # skip header
        jreader.next()
        # count rows
        for rowcount, _ in enumerate(jreader):
            pass
    return rowcount

def main():
    usage = 'usage: %prog [opt] swissprot_annotated_results '\
        'miso_summary output_file'\
        '\nThree arguments must be specified in command line:\n'\
        '1) swissprot_annotated_results (inclusion exons) '\
        '(output from annotate_genes_with_swissprot.py)\n'\
        '2) swissprot_annotaetd_results (exclusion exons) '\
        '(output from annotate_genes_with_swissprot.py)\n'\
        '3) miso_summary_file (bayes or T-Test).\n'\
        '4) inclusion protein summary (create_dna_protein_summary_file.py '\
            'output).\n'\
        '5) exclusion protein summary.(create_dna_protein_summary_file.py '\
            'output).\n'\
        '6) Output file\n'
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
    parser.add_option('--ylabel', dest='ylabel',
                      default='Counts',
                      help='ylabel of plot, defaults to "Counts"')
    parser.add_option('--title', dest='title',
                      default='Summary of UnitProt Features',
                      help='Title of plot. \n'\
                        'Default "Summary of UnitProt Features')
    parser.add_option('--normalize', dest='normalize',
                      default=True,
                      help='Normalize values by number of exons? True/False. '\
                        'Default True.')
    parser.add_option('--label1', dest='label1',
                      default='Inclusion',
                      help='Legend label 1: default "Inclusion"')
    parser.add_option('--label2', dest='label2',
                      default='Exclusion',
                      help='Legend label 2: default "Exclusion"')
    (options, args) = parser.parse_args()
    if len(args) != 6:
        print args
        print 'Incorrect number of args specified.\n%s' %usage
        sys.exit()
    # Define positional args
    swissprot_annot_fpath_incl = args[0]
    swissprot_annot_fpath_excl = args[1]
    miso_summary_fpath = args[2]
    incl_proteinsummary = args[3]    # for normalization
    excl_proteinsummary = args[4]    # for normalization
    output_fpath = args[5]
    # Parse options
    miso_event_colname = options.miso_event_colname
    samp1_colname = options.samp1_psi_colname
    samp2_colname = options.samp2_psi_colname
    summary_event_colname = options.summary_event_colname
    feature_colname = options.feature_colname
    ylabel = options.ylabel
    title = options.title
    label1 = options.label1
    label2 = options.label2
    # parse normalize, convert string to boolean
    normalize = options.normalize
    if normalize in ['True', 'TRUE', 'true', 'T', True]:
        normalize = True
    elif normalize in ['False', 'FALSE', 'false', 'F', False]:
        normalize = False
    else:
        print '%s must be True or Valse.' %normalize
    
    # Index miso summary, {misoevent: inclusion/exclusion}
    # allows separation of swissprot_annot_fpath_incl into incl or excl.
    inclexcl_dic = get_incl_excl_miso_events(miso_summary_fpath, 
                                             miso_event_colname, 
                                             samp1_colname, samp2_colname)
    
    # Determine number of exons by counting rows in protein summary
    # this is important if you want to normalize results.
    n_incl = get_exon_count(incl_proteinsummary)
    n_excl = get_exon_count(excl_proteinsummary)
    
    # Add n_incl and n_excl to label1 and label2
    label1 = ''.join([label1, ' (total DS events: %s)' %n_incl])
    label2 = ''.join([label2, ' (total DS events: %s)' %n_excl])
    
    incl_str, excl_str = get_features_dic_keys()
    # domains_dic {inclusion: {domain: }, exclusion: {domain: }}
    features_dic = {incl_str : {}, 
                    excl_str : {}}
    for incl_excl, swissprot_path in zip([incl_str, excl_str], 
                                         [swissprot_annot_fpath_incl, 
                                          swissprot_annot_fpath_excl]):
        features_subdic = get_features_dic(swissprot_path, 
                                           inclexcl_dic, 
                                           summary_event_colname, 
                                           feature_colname)
        features_dic[incl_excl].update(features_subdic[incl_excl])
        
    # Prepare features dic for bar plots
    ignore_list = ['VAR_SEQ', 'CHAIN']
    features_list, incl_count, excl_count = \
        prepare_dic_for_barplot(features_dic, ignore_list)
    
    if normalize:
        # normalize by incl excl counts
        incl_count = [100 * float(i) / n_incl for i in incl_count]
        excl_count = [100 * float(i) / n_excl for i in excl_count]
    
    # Plot bar plots
    plot_utils.plot_bar_plot(incl_count, excl_count, 
                             features_list, ylabel, title, 
                             label1, label2, width=0.4)
    
if __name__ == '__main__':
    main()