'''
Created on 2013-09-11

@author: jyeung
'''

import sys
import os
from multiprocessing import Queue
from t_test_miso_output import *
from group_miso_utils import *

if __name__ == '__main__':
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
    output_fname = sys.argv[5]
    
    # Define constants
    summary_fullpath = os.path.join(output_dir, output_fname)
    
    # Get sample names from textfile.
    group_1_samples = get_sample_names_from_file(group_1_samplenames_file)
    group_2_samples = get_sample_names_from_file(group_2_samplenames_file)
    
    # Create list of chromosomes.
    chr_list = create_chromo_list(prefix='chr')
    # chr_list = ['chr11']
    
    # Subset list for only those that contain miso outputs.
    group_1_samples = check_if_empty_dir(main_dir, group_1_samples, chr_list)
    group_2_samples = check_if_empty_dir(main_dir, group_2_samples, chr_list)
    
    # Init fnames dic
    fnames_dic = {}
    '''
    q = Queue()
    for chromo in chr_list:
        t_test_and_pickle(fnames_dic, chromo, output_dir, group_1_samples, group_2_samples, main_dir, q)
    '''
    
    chromo = 'chr5'
    
    q = Queue()
    t_test_and_pickle(fnames_dic, chromo, output_dir, group_1_samples, group_2_samples, main_dir, q, 10)
    fnames_dic.update(q.get())
    
    # Run on multiple threads.
    
    '''
    q = Queue()
    process_list = []
    for chromo in chr_list:
        print('Sending %s job to core...' %chromo)
        p = Process(target=t_test_and_pickle,
                    args=(fnames_dic, chromo, output_dir, 
                          group_1_samples, group_2_samples, 
                          main_dir, q))
        process_list.append(p)
        p.start()
    for chromo in chr_list:
        fnames_dic.update(q.get())
    
    # Wait for all threads to be done before continuing.
    for p in process_list:
        p.join()
    '''
        
    print('Completed %s jobs.' %len(chr_list))
    
    '''
    # Write fnames_dic as pickle file.
    pickle_filename = ''.join([output_fname, '_filenames_dic.pickle'])
    fnames_savepath = os.path.join(output_dir, pickle_filename)
    print('Saving filenames_dic.pickle to %s' %fnames_savepath)
    pickle_path = save_dic_as_pickle(fnames_dic, fnames_savepath)
    
    # Write information from pickle to textfile. 
    print('Writing information from pickle to textfile.')
    # Read pickle file to get fnames_dic
    fnames_dic = read_pickle(pickle_path)
    # Read and write to file. 
    read_pickle_write_to_file(summary_fullpath, chr_list, fnames_dic, 
                              filter_events=True)
    
    print('Summary file saved in: %s' %summary_fullpath)
    '''