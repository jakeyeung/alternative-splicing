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
from group_miso_utils import get_sample_names_from_file, create_chromo_list, \
    get_all_fnames, check_if_empty_dir, get_psi_dic_across_samples, \
    t_test_psi_info

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
    
    # Get sample names from textfile.
    group_1_samples = get_sample_names_from_file(group_1_samplenames_file)
    group_2_samples = get_sample_names_from_file(group_2_samplenames_file)
    all_samples = group_1_samples + group_2_samples
    
    # Create list of chromosomes.
    chr_list = create_chromo_list(prefix='chr')
    
    # Subset list for only those that contain miso outputs.
    group_1_samples = check_if_empty_dir(main_dir, group_1_samples, chr_list)
    group_2_samples = check_if_empty_dir(main_dir, group_2_samples, chr_list)
    
    # Redefine groups1 and groups2 to remove empty dirs.
    
    jchr = chr_list[0]    # For testing purposes.
    print jchr
    
    # Get list of AS events that need to be t-tested.
    group_1_fnames_list = get_all_fnames(group_1_samples, main_dir, jchr)
    group_2_fnames_list = get_all_fnames(group_2_samples, main_dir, jchr)
    master_fnames_list = group_1_fnames_list + group_2_fnames_list
    
    # Do t-test between the two groups. 
    for fname in master_fnames_list:
        # Get dictionary containing psi information for all samples.
        psi_info_dic, keynames = get_psi_dic_across_samples(fname, group_1_samples, 
                                                            group_2_samples, 
                                                            main_dir, jchr, 
                                                            output_dir)
        pval = t_test_psi_info(psi_info_dic)
        print pval
        raw_input()
        
if __name__ == '__main__':
    main()
    