'''
Created on 2014-01-21

@author: jyeung

Compare 3 cohorts, show overlap of
differentially expressed isoforms.
'''

import sys
import csv
from optparse import OptionParser
from utilities import plot_utils

def get_events(cohort_filename, colname='event_name'):
    '''
    Inputs:
        filename containing differentially expressed isoforms
        colname: colname containing miso ID. 
        colname defaults to None
    Outputs:
        list containing differentially expressed isoforms (genes?) 
    '''
    
    events = []
    with open(cohort_filename, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        header = myreader.next()
        for row in myreader:
            events.append(row[header.index(colname)])
    return set(events)

def parse_labels_csv(labels_option):
    '''
    Input:
        labels_option: 110,0.05,0.03;110,-0.01,0.03;...
    Output:
        labels_dic: {labelid: (x_adj, yadj), ...}
    
    Parse labels option (--adj_labels option flag)
    and return a dictionary of form:
    {labelid: (x_adj, yadj), ...}
    
    This dictionary is fed to plot function in order
    to adjust position of text labels.
    '''
    labels_dic = {}
    
    # parse labels option, first split by semicolon
    labels_ops_sc = labels_option.split(';')
    # expected form: labelid,x_adj,y-adj
    # split by comma, then store each element in dic
    for label_op_csv in labels_ops_sc:
        label_op = label_op_csv.split(',')
        labels_id = label_op[0]
        x_adj = float(label_op[1])
        y_adj = float(label_op[2])
        if labels_id not in labels_dic:
            labels_dic[labels_id] = (x_adj, y_adj)
        else:    # else complain
            print 'Error: Duplicate labels option: %s' %labels_id
            sys.exit()
    return labels_dic

def index_event_gsymbol(event_gsymbol_dic, 
                        cohort_fname, 
                        event_colname='event_name', 
                        gene_colname='gsymbol'):
    '''
    Update a dic (allows many fnames to be used)
    to get dic:
    {event: gsymbol, ...}
    '''
    with open(cohort_fname) as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        header = myreader.next()
        for row in myreader:
            event = row[header.index(event_colname)]
            gene = row[header.index(gene_colname)]
            if event not in event_gsymbol_dic:
                event_gsymbol_dic[event] = gene
    return event_gsymbol_dic

def main():
    usage = 'usage: %prog [opts] cohort1 cohort2 cohort3\n'\
        'Two args must be specified in commandline: \n'\
        '1,2,3) Text file containing differentially expressed isoforms.\n'\
        'Press -h or --help for option parameter information.'
    parser = OptionParser(usage=usage)
    parser.add_option('-1', '--colname1', dest='events_colname1',
                      default='event_name',
                      help='Colname containing events in cohort1')
    parser.add_option('-2', '--colname2', dest='events_colname2',
                      default='event',
                      help='Colname containing events in cohort2')
    parser.add_option('-3', '--colname3', dest='events_colname3',
                      default='event',
                      help='Colname containing events in cohort3')
    parser.add_option('--label1', dest='label1',
                      default='label1',
                      help='Venn diagram label for cohort1')
    parser.add_option('--label2', dest='label2',
                      default='label2',
                      help='Venn diagram label for cohort2')
    parser.add_option('--label3', dest='label3',
                      default='label3',
                      help='Venn diagram label for cohort3')
    parser.add_option('-t', '--title', dest='title',
                      default='Plot title',
                      help='Plot title.')
    parser.add_option('-a', '--adj_labels', dest='adj_labels',
                      default=None,
                      help='Adj labels, use CSVs separated by semi colons.\n'\
                        'Example:011,0.03,-0.03;110,-0.05,0.05\n'\
                        'Adjusts 011 label by +0.03 in x, -0.03 in y\n'\
                        'and adjust 110 label by -0.05 in x, 0.05 in y.\n'\
                        'Adjusting is usually only required for 3 set venns')
    # Parse options
    (options, args) = parser.parse_args()
    
    if len(args) != 3 and len(args) != 2:
        print usage
        sys.exit()
        
    if len(args) == 3:
        three_cohorts = True
    else:
        three_cohorts = False
    
    cohort1_fname = args[0]
    cohort2_fname = args[1]
    if three_cohorts:
        cohort3_fname = args[2]
    
    # Define relevant colnames of cohorts from options
    colname1 = options.events_colname1
    label1 = options.label1
    colname2 = options.events_colname2
    label2 = options.label2
    if three_cohorts:
        colname3 = options.events_colname3
        label3 = options.label3
    title = options.title
    
    # convert adj_labels option to a dic to be fed into plot function
    if options.adj_labels is not None:
        adj_labels_dic = parse_labels_csv(options.adj_labels)
    else:
        adj_labels_dic = options.adj_labels
    
    # get set of events
    cohort1_events = get_events(cohort1_fname, colname=colname1)
    cohort2_events = get_events(cohort2_fname, colname=colname2)
    if three_cohorts:
        cohort3_events = get_events(cohort3_fname, colname=colname3)
    
    # store dic of events:gsymbol
    event_gsymbol_dic = {}
    event_gsymbol_dic = index_event_gsymbol(event_gsymbol_dic, 
                                            cohort1_fname, event_colname=colname1, 
                                            gene_colname='gsymbol')
    
    if three_cohorts:
        print 'Overlapping genes:'
        for g in (cohort1_events & cohort2_events & cohort3_events):
            print event_gsymbol_dic[g]
        plot_utils.plot_three_set_venn(cohort1_events, 
                                       cohort2_events, 
                                       cohort3_events,
                                       adj_labels_dic, 
                                       mycolors=('c', 'y', 'g'), 
                                       mylabels=[label1, 
                                                 label2, 
                                                 label3],
                                       title=title)
    
    else:
        print 'Overlapping genes:'
        for g in (cohort1_events & cohort2_events):
            print event_gsymbol_dic[g]
        plot_utils.plot_two_set_venn(cohort1_events, 
                                     cohort2_events, 
                                     mycolors=('c', 'y'),
                                     mylabels=[label1, 
                                               label2],
                                     title=title)
    
if __name__ == '__main__':
    main()