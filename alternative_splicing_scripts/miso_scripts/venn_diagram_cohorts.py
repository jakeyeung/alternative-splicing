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
    # Parse options
    (options, args) = parser.parse_args()
    
    if len(args) != 2:
        print usage
        sys.exit()
    
    cohort1_fname = args[0]
    cohort2_fname = args[1]
    
    # Define relevant colnames of cohorts from options
    colname1 = options.events_colname1
    colname2 = options.events_colname2
    
    # get set of events
    cohort1_events = get_events(cohort1_fname, colname=colname1)
    cohort2_events = get_events(cohort2_fname, colname=colname2)
    
    print cohort1_events & cohort2_events
    
    plot_utils.plot_two_set_venn(cohort1_events, cohort2_events, 
                                 mylabels=['Xenografts: 331, 331R', 
                                           'Patient-Tumour Samples'])

if __name__ == '__main__':
    main()