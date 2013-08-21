'''
Created on 2013-07-23

@author: jyeung

Some functions to help out group_miso_results.py
'''

import os
import csv
import re
from scipy import stats
import numpy as np

def get_info_from_miso(psi_median_str, sample_name_str, psi_info_dic, 
                       main_dir, samp, chromo, fname):
    '''
    Checks if file exists, if exists, then add info to psi_info_dic
    '''
    # Create filepath
    file_path = os.path.join(main_dir, samp, chromo, fname)
    if os.path.exists(file_path):
        with open(file_path, 'rb') as readfile:
            reader = csv.reader(readfile, delimiter='\t')
            # Skip first two rows
            reader.next()
            '''
            TODO: Get header information from second row.
            '''
            reader.next()
            psi_value_list = []
            for row in reader:
                # First column (row[0]) is psi value,
                # Second column (row[1]) is log_score.
                # Psi value is comma separated, split it when take
                # the first value, which is the first isoform (inclusion ratio)
                psi_value_list.append(float(row[0].split(',')[0]))
            # keynames[0] should be psi_median_str from t_test_as_events()
            psi_info_dic[psi_median_str].append(np.median(psi_value_list))
            psi_info_dic[sample_name_str].append(samp)
    print(psi_info_dic)
    raw_input()
    return psi_info_dic
                
    
def t_test_as_events(master_fnames_list, group_1_samplenames, 
                      group_2_samplenames, main_dir, chromo, output_dir):
    '''
    Look at each as event across all samples between the two groups, then
    do a t-test to see if they are indeed different.
    
    Write this information to file. 
    '''
    # Define keyname constants
    psi_median_str = 'psi_median'
    sample_name_str = 'sample_name'
    counts_00_str = 'counts_00'
    counts_10_str = 'counts_10'
    counts_01_str = 'counts_01'
    counts_11_str = 'counts_11'
    assigned_counts_0_str = 'assigned_counts_0'
    assigned_counts_1_str = 'assigned_counts_1'
    keynames = [psi_median_str, sample_name_str, counts_00_str, counts_10_str, 
                counts_01_str, counts_11_str, 
                assigned_counts_0_str, assigned_counts_1_str]
    
    for fname in master_fnames_list:
        # Get psi information from each group as dictionary
        psi_info_dic = {}
        # Initialize keynames with empty list.
        for k in keynames:
            psi_info_dic[k] = []
        
        for samp1 in group_1_samplenames:
            file_dir = os.path.join(main_dir, samp1, chromo)
            # Get info from file
            psi_info_dic = get_info_from_miso(psi_median_str, sample_name_str, 
                                              psi_info_dic, 
                                              main_dir, samp1, 
                                              chromo, fname)
    print psi_info_dic
            
    
    
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