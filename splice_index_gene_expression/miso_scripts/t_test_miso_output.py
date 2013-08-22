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
from group_miso_utils import get_sample_names_from_file, create_chromo_list, \
    get_all_fnames, check_if_empty_dir, get_psi_dic_across_samples, \
    t_test_psi_info, save_dic_as_pickle, make_dir

def read_pickle_write_to_file(summary_fullpath, chr_list, fnames_dic, output_dir):
    '''
    Open a summary textfile, then individually open a pickle and write the 
    contents to file. 
    '''
    with open(summary_fullpath, 'wb') as writefile:
        master_fnames_list = fnames_dic[chr_list[0]]
        for fname in master_fnames_list:
            pickle_path = os.path.join(output_dir, chr_list[0], fname)
            pass
        
def t_test_and_pickle(fnames_dic, chromo, output_dir, group_1_samples, group_2_samples, 
                      main_dir):
    '''
    Combines several modules together into one so that the process
    can be easily multithreaded. 
    
    Return a dictionary containing chromosomes as keynames as fnames as values. 
    '''
    # Define constants
    pval_str = 'pval'
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
    
    # Do t-test between the two groups. 
    fnames_pickled_list = []
    for fname in master_fnames_list:
        # Get dictionary containing psi information for all samples.
        psi_info_dic, _ = get_psi_dic_across_samples(fname, 
                                                     group_1_samples, 
                                                     group_2_samples, 
                                                     main_dir, chromo, 
                                                     output_dir)
        psi_info_dic[pval_str] = t_test_psi_info(psi_info_dic)
        # Save dictionary as a pickle file.
        fnames_pickled_list.append(save_dic_as_pickle(psi_info_dic, chromo, 
                                                      fname, output_dir))
    # save fnames list to output dic
    if chromo not in fnames_dic:
        fnames_dic[chromo] = fnames_pickled_list
    else:
        print('Warning, overwriting fnames_list in %s' %chromo)
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
    
    # Select a chromo for testing
    jchr = chr_list[0]    # For testing purposes.
    print jchr
    
    fnames_dic = t_test_and_pickle(fnames_dic, chr_list, output_dir, 
                                   group_1_samples, group_2_samples, main_dir)
    print('T-tested all events in %s' %jchr)
    
    print('Writing information from pickle to textfile.')
    read_pickle_write_to_file(summary_fullpath, jchr, fnames_dic, output_dir)
        
        
if __name__ == '__main__':
    main()
    