'''
Created on 2013-08-21

@author: jyeung

Given a list of sample names between two classes (NEPC and PC, for example),
do a t-test between the two classes for all alternative splice events
(sample size may therefore vary) to test whether the PSI values differ
between the two classes. 

Keep in mind:
- not all AS events may be detected
- record log_score and confidence interval of each sample of each event 
    (miso gives probability distributions of psi rather than one value).
- re-use functions from group_miso_utils.py and group_miso_results.py
- use parallel processing for speed.
    -> naive strategy: create t-test file for each event in parallel,
        then consolidate all t-test files into one summary file?
    -> parallelize by chromosome. 
'''

import sys
import os
import csv
from multiprocessing import Process, Queue
from group_miso_utils import get_sample_names_from_file, create_chromo_list, \
    get_all_fnames, check_if_empty_dir, get_psi_dic_across_samples, \
    t_test_psi_info, save_dic_as_pickle, make_dir, read_pickle, get_psi_dic_keynames

def read_pickle_write_to_file(summary_fullpath, chr_list, fnames_dic, output_dir):
    '''
    Open a summary textfile, then individually open a pickle and write the 
    contents to file. 
    '''
    # Get keynames found in pickle file.
    # Each keyname will be a row written to file.
    _, psi_median_str, log_score_str, sample_name_str, \
        counts_00_str, counts_10_str, counts_01_str, counts_11_str, \
        assigned_counts_0_str, assigned_counts_1_str, \
        percent_accepted_str, group_str, pval_str, event_str \
            = get_psi_dic_keynames(full_keynames=True)
    
    writecount = 0
    with open(summary_fullpath, 'wb') as writefile:
        writer = csv.writer(writefile, delimiter='\t')
        # Write header
        header = [event_str, pval_str, sample_name_str, group_str, 
                         counts_00_str, counts_10_str, counts_01_str, 
                         counts_11_str, assigned_counts_0_str, 
                         assigned_counts_1_str, psi_median_str, percent_accepted_str, 
                         log_score_str]
        writer.writerow(header)
        
        for chromo in chr_list:
            pickle_fullpath_list = fnames_dic[chromo]
            for pickle_path in pickle_fullpath_list:
                psi_info_dic = read_pickle(pickle_path)
                row = []
                for key in header:
                    '''
                    # Dic contains both lists and strings.
                    But we want to only have one column per
                    keyvalue. Therefore, we collapse lists 
                    into comma separated values (CSV).
                    '''
                    if isinstance(psi_info_dic[key], basestring):
                        row.append(psi_info_dic[key])
                    else:
                        # It is not a string, assumes it is a list
                        # then join by comma.
                        row.append(','.join(psi_info_dic[key]))
                    writer.writerow(row)
                    writecount += 1
    return writecount
                
        
def t_test_and_pickle(fnames_dic, chromo, output_dir, group_1_samples, group_2_samples, 
                      main_dir):
    '''
    Combines several modules together into one so that the process
    can be easily multithreaded. 
    
    Return a dictionary containing chromosomes as keynames as fnames as values. 
    '''
    # Define constants
    pval_str = 'pval'
    event_str = 'event'
    # Define output dic
    fnames_dic = {}
    
    # Create directory to store pickled dictionary.
    make_dir(os.path.join(output_dir, chromo))
    
    '''
    # Get list of AS events that need to be t-tested.
    # Run the function on the lists separately to ensure
    # that each list contains at least one element.
    # This means our master_fnames_list is guaranteed to
    # have one sample in each group. 
    '''
    group_1_fnames_list = get_all_fnames(group_1_samples, main_dir, chromo)
    group_2_fnames_list = get_all_fnames(group_2_samples, main_dir, chromo)
    master_fnames_list = group_1_fnames_list + group_2_fnames_list
    # master_fnames_size = len(master_fnames_list)
    # Do t-test between the two groups. 
    fnames_pickled_list = []
    count = 0
    for fname in master_fnames_list:
        count += 1
        # Get dictionary containing psi information for all samples.
        psi_info_dic, _ = get_psi_dic_across_samples(fname, 
                                                     group_1_samples, 
                                                     group_2_samples, 
                                                     main_dir, chromo, 
                                                     output_dir)
        # Add pval and event to dic
        psi_info_dic[pval_str] = t_test_psi_info(psi_info_dic)
        # Remove .miso from fname to get event name. 
        psi_info_dic[event_str] = fname.split('.')[0]    
        # Save dictionary as a pickle file.
        # add .pickle to fname
        pickled_fname = ''.join([fname, '.pickle'])
        output_fullpath = os.path.join(output_dir, chromo, pickled_fname)
        fnames_pickled_list.append(save_dic_as_pickle(psi_info_dic, 
                                                      output_fullpath))
    # save fnames list to output dic
    if chromo not in fnames_dic:
        fnames_dic[chromo] = fnames_pickled_list
    else:
        print('Warning, overwriting fnames_list in %s' %chromo)
    print('T-tested all events in %s' %chromo)
    return fnames_dic
    
def main():
    '''
    Requires sample_list textfile, miso directory (where psi values are),
    and output_directory (where t-test results will be stored) 
    '''
    if len(sys.argv) < 3:
        print('Requires textfile of samplenames, miso_output directory,'\
              'and output_directory to be specified in command line.')
        sys.exit()
    group_1_samplenames_file = sys.argv[1]
    group_2_samplenames_file = sys.argv[2]
    main_dir = sys.argv[3]    # Where miso outputs results
    output_dir = sys.argv[4]
    
    # Define constants
    summary_fullpath = os.path.join(output_dir, 'summary_of_t_tests.txt')
    
    # Get sample names from textfile.
    group_1_samples = get_sample_names_from_file(group_1_samplenames_file)
    group_2_samples = get_sample_names_from_file(group_2_samplenames_file)
    
    # Create list of chromosomes.
    chr_list = create_chromo_list(prefix='chr')
    
    # Subset list for only those that contain miso outputs.
    group_1_samples = check_if_empty_dir(main_dir, group_1_samples, chr_list)
    group_2_samples = check_if_empty_dir(main_dir, group_2_samples, chr_list)
    
    # Init fnames dic
    fnames_dic = {}
    
    # Run on multiple threads.
    result_queue = Queue()
    for chromo in chr_list:
        print('Sending %s job to core...' %chromo)
        Process(target=t_test_and_pickle,
                args=(fnames_dic, chromo, output_dir, 
                      group_1_samples, group_2_samples, 
                      main_dir)).start()
    # Wait for all threads to be done before continuing.
    result_queue.get()    
    print('Completed %s jobs.' %len(chr_list))
    '''
    fnames_dic = t_test_and_pickle(fnames_dic, jchr, output_dir, 
                                   group_1_samples, group_2_samples, main_dir)
    '''
    # Write fnames_dic as pickle file.
    fnames_savepath = os.path.join(output_dir, 'filenames_dic')
    pickle_path = save_dic_as_pickle(fnames_dic, fnames_savepath)
    
    # Write information from pickle to textfile. 
    print('Writing information from pickle to textfile.')
    # Read pickle file to get fnames_dic
    fnames_dic = read_pickle(pickle_path)
    # Read and write to file. 
    read_pickle_write_to_file(summary_fullpath, chr_list, fnames_dic, output_dir)
        
        
if __name__ == '__main__':
    main()
    