'''
Created on 2013-07-23

@author: jyeung

Some functions to help out group_miso_results.py
'''

import os
import csv
import re

def make_dir(path):
    '''
    Make a directory if it doesn't exist.
    '''
    if not os.path.exists(path):
        os.makedirs(path)
    return path
    
def write_combined_miso_header(sample_dir_names_list, main_dir, 
                               chromo, miso_output, csv_write_obj):
    '''
    For given miso output file (.miso), find that .miso file from a list of
    samples, then write a common header by keeping everything the same
    but updating the percent_accept, count, and assigned_count.
    '''
    # Define constants found in header.
    percent_accept_index = 5
    counts_index = 7
    assigned_counts_index = 8
    # Initialize percent_accept, counts and assigned_counts
    percent_accept_list = []
    counts_00_list = []
    counts_10_list = []
    counts_01_list = []
    counts_11_list = []
    assigned_counts_0_list = []
    assigned_counts_1_list = []
    for samp_dir in sample_dir_names_list:
        with \
            open(os.path.join(main_dir, samp_dir, 
                              chromo, miso_output), 'rb') as readfile:
            
            reader = csv.reader(readfile, delimiter='\t')
            header = reader.next()
            # Get value of certain elements in header. 
            percent_str = header[percent_accept_index]
            # Get only the percent by taking the float from right
            # of equal sign.
            percent_accept_list.append(float(percent_str.split('=')[1]))
            
            # get counts and use regex to extract the counts:
            counts_str = header[counts_index]
            assigned_str = header[assigned_counts_index]
            # Find (0,0):, (1,0):, (0,1):, (1,1): and get integers after it:
            for jpat, jlist, jstr, jtype in zip(['0,0', '1,0', 
                                                 '0,1', '1,1',
                                                 '0', '1'], 
                                                [counts_00_list, 
                                                 counts_10_list, 
                                                 counts_01_list, 
                                                 counts_11_list,
                                                 assigned_counts_0_list,
                                                 assigned_counts_1_list],
                                                [counts_str]*4 + [assigned_str]*2,
                                                ['c']*4 + ['ass_c']*2):
                # Create regex expression to search
                if jtype == 'c':
                    reg_exprs = ''.join(['(?<=\(', jpat, '\)\:)\d+'])
                else:
                    reg_exprs = ''.join(['(?<=', jpat, '\:)\d+'])
                match = re.search(reg_exprs, jstr)
                # match = re.search('(?<=\(0,0\)\:)\d+', counts_str)
                if match:
                    try:
                        jlist.append(int(match.group(0)))
                    except ValueError:
                        print('Error: could not convert %s to int.' %match.group(0))
                else:
                    jlist.append(0)
    # Modify the last saved header with percent_accept, 
    # counts, assigned_counts and write as first line.
    # Modify percent_accept, counts and assigned_counts:
    header[percent_accept_index] = 'percent_accept=%.2f' \
        %(float(sum(percent_accept_list)) / len(percent_accept_list))
    header[counts_index] = 'counts=(0,0):%s,(1,0):%s,(0,1):%s,(1,1):%s' \
        %(sum(counts_00_list), sum(counts_10_list), 
          sum(counts_01_list), sum(counts_11_list))
    header[assigned_counts_index] = 'assigned_counts=0:%s,1:%s' \
        %(sum(assigned_counts_0_list), sum(assigned_counts_1_list))
    csv_write_obj.writerow(header)
    return None

def write_combined_psi_logscore(sample_dir_names_list, main_dir, 
                               chromo, miso_output, csv_write_obj):
    '''
    For a given miso_output file (that should be common across the
    samples of interest), find the sampled_psi and log_score for
    all samples in sample_dir_names_list. Write this to csv_write_obj.
    '''
    # Define column headers.
    sampledpsi_str = 'sampled_psi'    # First column is sampled_psi
    logscore_str = 'log_score'    # Second column is log_score
    # Write column headers to write_obj.
    csv_write_obj.writerow([sampledpsi_str, logscore_str])
    sampled_psi_count = 0
    for samp_dir in sample_dir_names_list:
        with \
            open(os.path.join(main_dir, samp_dir, 
                              chromo, miso_output), 'rb') as readfile:
            reader = csv.reader(readfile, delimiter='\t')
            # Skip first row and second row.
            reader.next()
            reader.next()
            # Itereate rows containing the sampled_psi and log_score values.
            for row in reader:
                csv_write_obj.writerow(row)
                sampled_psi_count += 1
    return sampled_psi_count