'''
Created on 2013-09-09

@author: jyeung

Strategy:
Load a set of putative switch-AS-events as dic containing expected psi val 
    for NEPC/PC
Read line by line in the miso_sample and for ones that match putative 
switch-AS-events:
    -Record PSI value for that sample
    -Check which side of the switch-AS-event (PC or NEPC) then
    -make a prediction of whether it should be PC or NEPC? 
'''


import sys
import csv
import re


class events_read_obj(object):
    '''
    For reading known as events.
    '''
    def __init__(self, path, has_header=True):
        self.path = path
        if has_header:
            self.has_header = True
        else:
            self.has_header = False
    
    def __enter__(self):
        self.readfile = open(self.path, 'rb')
        self.reader = csv.reader(self.readfile, delimiter='\t')
        if self.has_header:
            self.header = self.reader.next()
        else:
            self.header = None
    
    def __exit__(self, exittype, exitvalue, exittraceback):
        self.readfile.close()
        
def get_col_index(jstr, header):
    '''
    Get column index, given a header and column name.
    '''
    myindex = header.index(jstr)
    return(myindex)
        
def convert_csv_to_list(mycsv, convert_to_int=True):
    '''
    Given csv values, splits by comma and returns list.
    If convert_to_int, then will try to convert values to int.
    '''
    # Initialize constants
    mylist = mycsv.split(',')
    if convert_to_int:
        mylist = [int(i) for i in mylist]
    return mylist

def above_avg_min_counts(row, header, avg_min_counts_thres=1):
    '''
    Check if the row in events file contains counts_01 greater
    than min_counts by calculating average min_counts.
    '''
    counts01_index = get_col_index('counts_01', header)
    counts01_list = convert_csv_to_list(row[counts01_index], 
                                        convert_to_int=True)
    # Calculate average min counts
    my_avg_min_count = float(sum(counts01_list)) / len(counts01_list)
    if my_avg_min_count < avg_min_counts_thres:
        above_threshold = False
    else:
        above_threshold = True
    return above_threshold
    
def above_min_psi_diff(row, header, psi_diff_thres=0.4):
    '''
    Check if column abs_delta_psi is greater than threshold or not.
    '''
    abs_delta_psi_index = get_col_index('abs_delta_psi', header)
    abs_delta_psi = float(row[abs_delta_psi_index])
    if abs_delta_psi < psi_diff_thres:
        above_threshold = False
    else:
        above_threshold = True
    return above_threshold

def get_header_strings():
    '''
    Get some predefined colnames used as 
    sub-key values.
    '''
    # Define constants
    gsymbol_str = 'gsymbol'
    sample_name_str = 'sample_name'
    group_str = 'group'
    counts01_str = 'counts_01'
    psi_medians_str = 'psi_median'
    delta_psi_str = 'delta_psi'
    header_list = [gsymbol_str, sample_name_str, group_str, 
                   counts01_str, psi_medians_str, delta_psi_str]
    return header_list

def get_output_header():
    '''
    Run get_header_strings, then add event, sample_psi_mean and sample_counts
    '''
    # Define str constants
    sample_psi_mean_str = 'sample_psi_mean'
    sample_counts_str = 'sample_counts'
    
    prelim_header = get_header_strings()
    header = prelim_header + [sample_psi_mean_str, sample_counts_str]
    
    return header

def init_key_to_dic(mykey, mydic):
    '''
    Checks if key exists, then adds key to dic and creates
    an empty dictionary to initialize for further writing.
    '''
    if mykey not in mydic:
        mydic[mykey] = {}
    else:
        print('Warning, %s exists in dic, ovewriting...' %mykey)
    return(mydic)

def write_info_to_dic(row, event, header, colnames, dic):
    '''
    Given an empty dic, add information in the form of
    key:values (colname:row[col_index])
    '''
    # Initializek key
    dic = init_key_to_dic(event, dic)
    # Add values to initialized key dic
    for c in colnames:
        dic[event][c] = row[header.index(c)]
    return dic

def write_event_to_dic(row, header, dic):
    '''
    Write necessary info to dic such as
    miso-event to key, gsymbol, groups, psi_medians, delta_psi, counts01
    '''
    my_colnames = get_header_strings()
    events_index = get_col_index('event', header)
    event = row[events_index]
    dic = write_info_to_dic(row, event, header, my_colnames, dic)
    return dic

def iterate_filter_events(events_obj):
    '''
    Iterate rows, retaining only rows that pass filter.
    '''
    events_dic = {}
    rowcount = 0
    for row in events_obj.reader:
        if above_avg_min_counts(row, events_obj.header, avg_min_counts_thres=1) and \
            above_min_psi_diff(row, events_obj.header, psi_diff_thres=0.4):
            write_event_to_dic(row, events_obj.header, events_dic)
            rowcount += 1
    print('%s events indexed into dictionary.' %rowcount)
    return events_dic

def index_high_conf_as_events(as_events_path, min_psi_diff=0.4):
    '''
    Look through as_events file of top switch-AS events,
    filter out low counts_01 events and filtero ut ones below
    min_psi_diff
    '''
    events_obj = events_read_obj(as_events_path)
    with(events_obj):
        indexed_events = iterate_filter_events(events_obj)
    return indexed_events

def get_sample_info(miso_summary_path):
    '''
    Read miso_summary file, return a list of tuples containing
    (eventname, miso_posterior_mean, counts)
    '''
    # Define constants
    event_name_index = 0
    psi_mean_index = 1
    counts_index = 5
    
    # store as info tuple list
    tup_list = []
    
    with open(miso_summary_path, 'rb') as misofile:
        misoreader = csv.reader(misofile, delimiter='\t')
        misoreader.next()    # Assume has header, skip it.
        for row in misoreader:
            tup_list.append((row[event_name_index], row[psi_mean_index], row[counts_index]))
    return tup_list

def append_element_to_index_dic(index_dic, key, subkey, subval):
    '''
    Append information to a subdic
    Assumes a dic within a dic.
    '''
    try:
        index_dic[key][subkey] = subval
    except KeyError:
        pass
    return index_dic

def append_sample_info_to_index_dic(sample_info, index_dic):
    '''
    Add sample_info to dictionary, then we can write it
    to file. 
    '''
    for tup in sample_info:
        samp_event = tup[0]
        samp_psi_mean = tup[1]
        samp_counts = tup[2]
        for info_str, info in \
            zip(['sample_psi_mean', 'sample_counts'], 
                [samp_psi_mean, samp_counts]):
            index_dic = \
                append_element_to_index_dic(index_dic, samp_event, 
                                         info_str, info)
    return index_dic

def write_header_to_file(header, write_obj):
    '''
    Get a list of strings called header, write it as first line to write_obj.
    '''
    write_obj.writerow(['event'] + header)
    
def above_counts_threshold(counts, threshold=3):
    '''
    Is counts (0,1) or (1,0) above a certain threshold?
    '''
    try:
        counts01 = int(re.search('(?<=\(0,1\)\:)\w+', counts).group(0))
    except AttributeError:
        counts01 = 0
    try:
        counts10 = int(re.search('(?<=\(1,0\)\:)\w+', counts).group(0))
    except AttributeError:
        counts10 = 0
    count_total = counts01 + counts10
    if count_total < threshold:
        above_counts = False
    else:
        above_counts = True
    
    return above_counts
    
def write_embedded_dic_to_file(index_dic, header, write_obj):
    '''
    Take embedded key with each key as row, each subkey as column.
    '''
    writecount = 0
    miss_count = 0
    for key, subdic in index_dic.iteritems():
        try:
            counts = subdic['sample_counts']
            no_samp = False
        except KeyError:
            no_samp = True
        if above_counts_threshold(counts, threshold=3) or no_samp==True:
            wrow = []
            wrow.append(key)    # Make event the first element in list.
            for h in header:
                try:
                    wrow.append(subdic[h])
                except KeyError:
                    miss_count += 1
            write_obj.writerow(wrow)
            writecount += 1
    print('%s events missed for some reason...' %miss_count)
    return writecount

def write_sample_info_and_dic_to_file(index_dic, sample_info, write_path):
    '''
    In write_path, write the index_dic used to find sample_info, then 
    for each event, write the info inside the index_dic as well as the
    sample_info for that event.
    '''
    with open(write_path, 'wb') as writepath:
        jwriter = csv.writer(writepath, delimiter='\t')
        # sample_info is tuple_list in form (event, psi_mean, counts)
        index_dic = append_sample_info_to_index_dic(sample_info, index_dic)
        # Write header
        out_header = get_output_header()
        write_header_to_file(out_header, jwriter)
        # Write embedded dic tofile.
        write_count = write_embedded_dic_to_file(index_dic, out_header, jwriter)
    print('%s rows written to %s.' %(write_count, write_path))
    return index_dic
                    
def main(as_events_path, miso_summary_path, write_path):
    as_events_dic = index_high_conf_as_events(as_events_path)
    tup_list = get_sample_info(miso_summary_path)
    appended_dic = write_sample_info_and_dic_to_file(as_events_dic, 
                                                     tup_list, 
                                                     write_path)

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('File containing AS events and file containing '\
              'miso_summary must be specified in command line.')
        sys.exit()
    as_events_path = sys.argv[1]
    miso_summary_path = sys.argv[2]
    write_path = sys.argv[3]
    main(as_events_path, miso_summary_path, write_path)