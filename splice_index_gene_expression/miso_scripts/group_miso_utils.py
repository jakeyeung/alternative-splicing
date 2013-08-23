'''
Created on 2013-07-23

@author: jyeung

Some functions to help out group_miso_results.py
'''

import sys
import os
import csv
import re
import pickle
from scipy import stats
import numpy as np

def read_pickle(pickle_path):
    with open(pickle_path, 'rb') as pkl_file:
        mydata = pickle.load(pkl_file)
    return mydata

def save_dic_as_pickle(dic, output_fullpath, protocol=-1):
    '''
    Saves dictionary to output path.
    
    Protocol -1 uses highest protocol available. 
    If protocol = 0, uses standard protocol 0.
    '''
    # Define constants
    pickle_str = '.pickle'
    
    # Add prefix 'pickle' to end of fullpath
    output_fullpath = ''.join([output_fullpath, pickle_str])
    
    with open(output_fullpath, 'wb') as output:
        pickle.dump(dic, output, -1)
    return output_fullpath

def split_psi_medians_into_two_lists(info_list, group_list):
    '''
    Given list of info (e.g. psi medians) and corresponding group list, 
    separate psi_medians into two separate lists. 
    '''
    # Check the lengths are the same for both lists
    if len(info_list) != len(group_list):
        print('info_list has length %s, not equal to length of '\
              'group_list, %s' %(len(info_list), len(group_list)))
        sys.exit()
    # Define possible values in group_list
    group_1_str = '1'    # means this is in group1
    group_2_str = '2'    # means this is in group2
    
    # Init output lists
    group_1_list = []
    group_2_list = []
    
    for i, g in enumerate(group_list):
        if g==group_1_str:
            group_1_list.append(info_list[i])
        elif g==group_2_str:
            group_2_list.append(info_list[i])
        else:
            print('Cannot tell if %s is in group_1 or group_2' %g)
            sys.exit()
    return group_1_list, group_2_list
    
def t_test_psi_info(psi_info_dic, psi_median_str='psi_median', 
                    group_str='group'):
    '''
    Given psi info dic which is a dictionary that 
    contains psi_median for all samples, do a t-test between
    group1 and group2 (group info found inside psi_info).
    
    Strings psi_median and group are used as keys to access
    information in dictionary.
    '''
    # Get lists from dic
    group_list = psi_info_dic[group_str]
    psi_median_list = psi_info_dic[psi_median_str]
    
    # Separate the psi_median_list into two groups
    psi_medians_group1, psi_medians_group2 = \
        split_psi_medians_into_two_lists(psi_median_list, group_list)
    
    # T-test the two groups
    if len(psi_medians_group1) > 1 and len(psi_medians_group2) > 1:
        _, pval = stats.ttest_ind(psi_medians_group1, psi_medians_group2)
    else:
        pval = 'NA'
    return pval
    
def get_percent_accepted_from_header(header, percent_accept_index = 5):
    '''
    Read the row, assumes it is the miso header which contains (by default):
        percent_accepted column at index 5.
        
    Output percent accept 
    '''
    # Get column names from index
    percent_str = header[percent_accept_index]
    # Get the right side of equal sign from e.g. "percent_accept=99"
    percent_accept = float(percent_str.split('=')[1])
    return percent_accept
    
def read_counts_from_miso_header(header, 
                                 counts_index = 7, 
                                 assigned_counts_index = 8):
    '''
    Read the row, assumes it is the miso header which contains (by default):
        counts (0,0), (1,0), (0,1), (1,1) at index 7
        assigned-counts (0, 1) at index 8.
        
    Output counts and assigned-counts to user. 
    '''
    # Get column names from index
    counts_str = header[counts_index]
    assigned_str = header[assigned_counts_index]
    
    '''
    # Get counts by finding (0,0):, (1,0):, (0,1):, (1,1)
    # and
    # Get assigned_counts by finding 0:, 1:
    '''
    # Find (0,0):, (1,0):, (0,1):, (1,1): and get integers after it:
    # Initialize variable
    counts_00 = []    # Reads not matching any isoform
    counts_10 = []    # Reads matching to inclusion isoform
    counts_01 = []    # Reads matching to exclusion isoform
    counts_11 = []    # Reads matching both isoforms (exon bodies)
    assigned_counts_0 = []    # Counts assigned to inclusion isoform.
    assigned_counts_1 = []    # Counts assigned to exclusion isoform.
    
    for jpat, jlist, jstr, jtype in \
        zip(['0,0', '1,0', '0,1', '1,1','0', '1'], 
            [counts_00, counts_10, counts_01, counts_11, assigned_counts_0, 
             assigned_counts_1],
            [counts_str]*4 + [assigned_str]*2,
            ['c']*4 + ['ass_c']*2):
        
        # Create regex expression to search
        if jtype == 'c':    # counts
            # Search digits after (0,0):
            reg_exprs = ''.join(['(?<=\(', jpat, '\)\:)\d+'])
        elif jtype == 'ass_c':    # assigned_counts, 'ass_c'
            # Search digits after 0:
            reg_exprs = ''.join(['(?<=', jpat, '\:)\d+'])
        else:
            print('Expected %s to be either "c"(counts) '\
                  'or "ass_c"(assigned_counts)')
            sys.exit()
        match = re.search(reg_exprs, jstr)
        # match = re.search('(?<=\(0,0\)\:)\d+', counts_str)
        if match:
            try:
                jlist.append(int(match.group(0)))
            except ValueError:
                ('Error: could not convert %s to int.' %match.group(0))
        else:
            jlist.append(0)
    '''
    Return first element from each jlist. I use lists instead of a variable
    so I could loop through each. 
    '''
    return counts_00[0], counts_11[0], counts_01[0], counts_11[0], \
            assigned_counts_0[0], assigned_counts_1[0]
            
def get_info_from_miso(psi_median_str, log_score_str, 
                       sample_name_str, 
                       counts_00_str, counts_10_str, 
                       counts_01_str, counts_11_str, 
                       assigned_counts_0_str, 
                       assigned_counts_1_str, 
                       percent_accepted_str,
                       group_str,
                       psi_info_dic, 
                       samp,
                       group_1_samplenames,
                       group_2_samplenames, 
                       file_path):
    '''
    Checks if file exists, if exists, then add info to psi_info_dic
    '''
    if os.path.exists(file_path):
        # Record whether sample is in group1 or group2.
        if samp in group_1_samplenames:
            psi_info_dic[group_str].append('1')
        elif samp in group_2_samplenames:
            psi_info_dic[group_str].append('2')
        else:
            print('Could not place %s in either group 1 or group 2.' %samp)
        
        with open(file_path, 'rb') as readfile:
            reader = csv.reader(readfile, delimiter='\t')
            '''
            First row contains percent_accepted and counts info.
            Second row is just column names (psi_value, log_score)
            
            So read first row as header, but skip second row.
            '''
            header = reader.next()
            reader.next()
            # Get percent_accepted values
            psi_info_dic[percent_accepted_str].append\
                (get_percent_accepted_from_header(header))
            # Get counts
            counts_00, counts_10, counts_01, counts_11, assigned_counts_0, \
                assigned_counts_1 = read_counts_from_miso_header(header)
            
            # Put counts into dic
            for key, var in zip([counts_00_str, counts_10_str, 
                                   counts_01_str, counts_11_str, 
                                   assigned_counts_0_str, 
                                   assigned_counts_1_str], 
                                  [counts_00, counts_10, counts_01, counts_11, 
                                   assigned_counts_0, assigned_counts_1]):
                psi_info_dic[key].append(var)
            
            psi_value_list = []
            log_score_list = []
            for row in reader:
                # First column (row[0]) is psi value,
                # Second column (row[1]) is log_score.
                # Psi value is comma separated, split it when take
                # the first value, which is the first isoform (inclusion ratio)
                psi_value_list.append(float(row[0].split(',')[0]))
                log_score_list.append(float(row[1]))
            # keynames[0] should be psi_median_str from t_test_as_events()
            psi_info_dic[psi_median_str].append(np.median(psi_value_list))
            psi_info_dic[log_score_str].append(np.median(log_score_list))
            psi_info_dic[sample_name_str].append(samp)
    else:    # File doesn't exist
        # print('%s does not exist for sample %s' %(file_path, samp))
        pass
    return psi_info_dic

def get_psi_dic_keynames(full_keynames=False):
    '''
    Gets keynames used in get_psi_dic_across_samples()
    '''
    # Define keyname constants
    psi_median_str = 'psi_median'
    log_score_str = 'log_score_median'
    sample_name_str = 'sample_name'
    counts_00_str = 'counts_00'
    counts_10_str = 'counts_10'
    counts_01_str = 'counts_01'
    counts_11_str = 'counts_11'
    assigned_counts_0_str = 'assigned_counts_0'
    assigned_counts_1_str = 'assigned_counts_1'
    percent_accepted_str = 'percent_accepted'
    group_str = 'group'
    # Create keynames from constants
    keynames = [psi_median_str, log_score_str, sample_name_str, counts_00_str, 
                counts_10_str, counts_01_str, counts_11_str, 
                assigned_counts_0_str, assigned_counts_1_str, 
                percent_accepted_str, group_str]
    if full_keynames==True:
        pval_str = 'pval'
        event_str = 'event'
        keynames.append(pval_str)
    else:
        pval_str = 'pval'
        event_str = 'event'
    
    return keynames, psi_median_str, log_score_str, sample_name_str, \
        counts_00_str, counts_10_str, counts_01_str, counts_11_str, \
        assigned_counts_0_str, assigned_counts_1_str, \
        percent_accepted_str, group_str, pval_str, event_str
    
def get_psi_dic_across_samples(fname, group_1_samplenames, 
                                group_2_samplenames, main_dir, chromo, 
                                output_dir):
    '''
    Look at each as event across all samples between the two groups, then
    do a t-test to see if they are indeed different.
    
    Write this information to file. 
    '''
    # Get keynames for psi_info_dic
    keynames, psi_median_str, log_score_str, sample_name_str, \
            counts_00_str, counts_10_str, counts_01_str, counts_11_str, \
            assigned_counts_0_str, assigned_counts_1_str, \
            percent_accepted_str, group_str, _, _ = get_psi_dic_keynames()
            
    # Get psi information from each group as dictionary
    psi_info_dic = {}
    # Initialize keynames with empty list.
    for k in keynames:
        psi_info_dic[k] = []
    for samp in group_1_samplenames + group_2_samplenames:
        file_dir = os.path.join(main_dir, samp, chromo, fname)
        # Get psi info from file
        psi_info_dic = get_info_from_miso(psi_median_str, log_score_str,
                                          sample_name_str, 
                                          counts_00_str, counts_10_str, 
                                          counts_01_str, counts_11_str, 
                                          assigned_counts_0_str, 
                                          assigned_counts_1_str, 
                                          percent_accepted_str,
                                          group_str,
                                          psi_info_dic, 
                                          samp, 
                                          group_1_samplenames,
                                          group_2_samplenames,
                                          file_dir)
    return psi_info_dic, keynames
    
def get_all_fnames(sample_dir_list, main_dir, chromo):
    '''
    For each sample,
    
    Looks into miso_dir/samp_dir/chromo and gets the filename of 
    every file in that directory. 
    
    Returns a list of all possible filenames across all samples.
    There will be no duplicates in this list. 
    '''
    master_fnames_list = []
    
    for samp_dir in sample_dir_list:
        file_dir = os.path.join(main_dir, samp_dir, chromo)
        fnames_list = os.listdir(file_dir)
        # Add files with ending .miso in list to master list.
        master_fnames_list += [f for f in fnames_list if f.endswith('.miso')]
        if len(fnames_list) == 0:
            print('Warning: no .miso files found in %s.' %samp_dir)
    # Remove dulpicates.
    master_fnames_list = list(set(master_fnames_list))
    return master_fnames_list
    
    
def create_chromo_list(prefix='chr'):
    '''
    Create list of chromosomes (chr1, chr2, chr3... chr22, chrX, chrY)
    Assumes prefix of chr unless otherwise specified.
    '''
    # Define constants
    chr_str = prefix
    # Create list of chromosome names corresponding to folders within sample dir
    chr_list = [''.join([chr_str, str(c)]) for c in range(1, 23) + ['X', 'Y']]
    return chr_list
    
def check_if_empty_dir(path, directory_list, suffix_dir_list):
    '''
    Given a path+directory+suffix_dir, check if there are any files
    in the directory. If no files, remove that directory
    from the list, return only directories that contain files.
    '''
    bad_dirs = []
    filter_count = 0
    for samp_dir in directory_list:
        # Check if chr1, chr2 folders even exist...
        mainpath = os.path.join(path, samp_dir)
        files_in_mainpath = os.listdir(mainpath)
        if len(files_in_mainpath) == 0:
            bad_dirs.append(samp_dir)
            filter_count += 1
        else:    # Dont go into loop unless there are actually folders
            for suf_dir in suffix_dir_list:
                fullpath = os.path.join(path, samp_dir, suf_dir)
                files_in_fullpath = os.listdir(fullpath)
                files_in_fullpath = \
                    [f for f in files_in_fullpath if f.endswith('.miso')]
                if len(files_in_fullpath) == 0:
                    bad_dirs.append(samp_dir)
                    filter_count += 1
                    break
    filtered_dir_list = [f for f in directory_list if f not in bad_dirs]
    print('%s total samples were filtered out because they were empty.' \
          %filter_count)
    return filtered_dir_list

def get_sample_names_from_file(sample_name_file):
    '''
    Open text file containing sample names on the columns.
    
    Extract first column containing sample names, while appending
    a prefix in front of each sample name.
    
    Return sample names as a list. 
    '''
    sample_list = []
    with open(sample_name_file, 'rb') as readfile:
        reader = csv.reader(readfile, delimiter='\t')
        for row in reader:
            sample_list.append(row[0])
    return sample_list

def get_larger_list_size(list_of_interest, current_size):
    '''
    Check if current length of output_fnames_list is 
    larger than largest_fnames_size, if it is larger
    then update largest_fnames_size
    '''
    if len(list_of_interest) > current_size:
        return len(list_of_interest)
    else:
        return current_size

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
                if jtype == 'c':    # counts
                    # Search digits after (0,0):
                    reg_exprs = ''.join(['(?<=\(', jpat, '\)\:)\d+'])
                else:    # assigned_counts
                    # Search digits after 0:
                    reg_exprs = ''.join(['(?<=', jpat, '\:)\d+'])
                match = re.search(reg_exprs, jstr)
                # match = re.search('(?<=\(0,0\)\:)\d+', counts_str)
                if match:
                    try:
                        jlist.append(int(match.group(0)))
                    except ValueError:
                            ('Error: could not convert %s to int.' %match.group(0))
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