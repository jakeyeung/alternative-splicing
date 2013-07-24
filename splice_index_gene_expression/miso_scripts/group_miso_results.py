'''
Created on 2013-07-23

@author: jyeung

MISO results contains numerous files and folders for each sample.
In our Mark Rubin cohort, we have 38 PCa samples and 4 NEPCa.
Let's create one summary directory containing all the 38 PCa samples and their events.
Secondly, we'll create similarly for the 4 NEPCa one summary directory.
'''


import sys
import csv
import os
from multiprocessing import Process
from group_miso_utils import \
    write_combined_miso_header, make_dir, write_combined_psi_logscore, \
    get_larger_list_size, get_sample_names_from_file, check_if_empty_dir


def open_miso_output_file(miso_outpath):
    '''
    TODO: THIS FUNCTION'S USELESS DELETE IT. 
    Open the raw output from miso pipeline and extract information from it.
    '''
    with open(miso_outpath, 'rb') as readfile:
        reader = csv.reader(readfile, delimiter='\t')
        print reader.next()
    
def get_intersection_of_files(main_dir, sample_dir_names_list, chromo):
    '''
    Itereate through list of sample_dir_names and obtain a set of .miso output
    files in miso_main_dir/sample_dir/chromo/ that are common through 
    all the sample directories that were iterated. 
    
    Also output total possible fnames by finding the list of largest length
    as we iterate through samples. This gives us an idea of how many get filtered out.
    '''
    # Initialize master_fnames_list
    largest_fnames_size = 0    # Total possible fnames.
    master_fnames_list = []
    samp_count = 0
    for samp_dir in sample_dir_names_list:
        chr_output_file_dir = os.path.join(main_dir, samp_dir, chromo)
        output_fnames_list = os.listdir(chr_output_file_dir)
        # Filter for .miso files only.
        output_fnames_list = [f for f in output_fnames_list if f.endswith('.miso')]
        if len(output_fnames_list) == 0:
            print('Warning: no .miso files found in %s.' %samp_dir)
        largest_fnames_size = get_larger_list_size(output_fnames_list, largest_fnames_size)
        '''
        If this is first iteration, make fnames_list the masters list.
        Otherwise, find intersection between masters list and fnames list
        '''
        if samp_count > 0:
            if len(master_fnames_list) > 0:
                master_fnames_list = \
                    list(set(master_fnames_list).intersection(output_fnames_list))
            elif len(master_fnames_list) == 0:
                # print('Warning: master_fnames_list is zero.')
                pass
            else:
                sys.exit('Cannot find length of master_fnames_list.')
        else:    # If this is first iteration, initialize master_fnames_list.
            master_fnames_list = output_fnames_list
        samp_count += 1
    # print('Found common %s events in %s.' %(len(master_fnames_list), chromo))
    return master_fnames_list, largest_fnames_size

def consolidate_miso_across_samples(main_dir, sample_dir_names_list, 
                                    chromo, master_fnames_list, 
                                    output_path):
    '''
    For each file in master_fnames_list:
    Read inside main_dir/sample_dir/chromo/misofile for all sample_dirs.
    1) read the header and write a combined header to the output_path
    (meaning adding up all the counts).
    2) write sampled_psi and log_score to output_path for each
    sample. 
    
    TODO: is it weird to put csv writer obj into a function?
    '''
    # This will be inefficient, because I will open each file twice...
    file_count = 0
    for f in master_fnames_list:
        # Construct fullpath for output file, make directory if does not exist.
        # Path will look like output_path/chr1/
        chr_out_path = make_dir(os.path.join(output_path, chromo))
        with open(os.path.join(chr_out_path, f), 'wb') as writefile:
            writer = csv.writer(writefile, delimiter='\t')
            # Write summary header of file (loops through all samples)
            write_combined_miso_header(sample_dir_names_list, main_dir, 
                                       chromo, f, writer)
            # Write sampled_psi and log_score for each sample:
            write_combined_psi_logscore(sample_dir_names_list,
                                        main_dir,
                                        chromo,
                                        f, writer)
            file_count += 1
    # print('%s files summarized for %s samples.' %(file_count, len(sample_dir_names_list)))

def get_fnames_consolidate_miso(main_dir, sample_dir_names_list, chromo, 
                                output_path):
    '''
    Simply combines the two functions get_intersection_of_files() and 
    consolidate_miso_across_samples(). This makes it easier to create 
    queues for multiprocessing.
    '''
    master_fnames_list, \
        total_possible = get_intersection_of_files(main_dir, 
                                                   sample_dir_names_list, 
                                                   chromo)
    # Create first header for each common file name.
    consolidate_miso_across_samples(main_dir, sample_dir_names_list,
                                    chromo, master_fnames_list,
                                    output_path)
    print('%s files (of %s total possible) consolidated for chromosome %s' \
          %(len(master_fnames_list), total_possible, chromo))
    return None

def main():
    main_dir = sys.argv[1]
    sample_dir_fullpath = sys.argv[2]
    output_path = sys.argv[3]
    
    # Define constants
    chr_str = 'chr'
    # Create list of chromosome names corresponding to folders within sample dir
    chr_list = [''.join([chr_str, str(c)]) for c in range(1, 23) + ['X', 'Y']]
    
    # Create list of sample directory names.
    sample_dir_names_list = get_sample_names_from_file(sample_dir_fullpath)
    
    # Subset list for only those that contain miso outputs.
    sample_dir_names_list = check_if_empty_dir(main_dir, sample_dir_names_list, chr_list)
    # sample_dir_names_list = sample_dir_names_csv.split(',')
    
    # Run on multiple threads.
    for chromo in chr_list:
        print('Sending %s job to core...' %chromo)
        Process(target=get_fnames_consolidate_miso,
                args=(main_dir, sample_dir_names_list, chromo, 
                      output_path)).start()

if __name__ == '__main__':
    main()