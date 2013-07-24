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
from group_miso_utils import write_headers_all_samples


def open_miso_output_file(miso_outpath):
    '''
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
    '''
    # Initialize master_fnames_list
    master_fnames_list = []
    samp_count = 0
    for samp_dir in sample_dir_names_list:
        print samp_dir
        chr_output_file_dir = os.path.join(main_dir, samp_dir, chromo)
        output_fnames_list = os.listdir(chr_output_file_dir)
        # Filter for .miso files only.
        output_fnames_list = [f for f in output_fnames_list if f.endswith('.miso')]
        if len(output_fnames_list) == 0:
            print('Warning: no .miso files found in %s.')
        '''
        If this is first iteration, make fnames_list the masters list.
        Otherwise, find intersection between masters list and fnames list
        '''
        if samp_count > 0:
            if len(master_fnames_list) > 0:
                master_fnames_list = \
                    list(set(master_fnames_list).intersection(output_fnames_list))
            elif len(master_fnames_list) == 0:
                print('Warning: master_fnames_list is zero.')
                pass
            else:
                sys.exit('Cannot find length of master_fnames_list.')
        else:    # If this is first iteration, initialize master_fnames_list.
            master_fnames_list = output_fnames_list
        samp_count += 1
    return master_fnames_list



def combine_miso_outputs_across_samples(main_dir, sample_dir_names_list, 
                                        chromo, master_fnames_list, 
                                        output_path):
    '''
    For each file in master_fnames_list:
    Read inside main_dir/sample_dir/chromo/misofile for all sample_dirs.
    1) read the header and write a combined header to the output_path
    (meaning adding up all the counts).
    2) write sampled_psi and log_score to output_path for each
    sample. 
    '''
    
    '''
    # Define constants found in header.
    percent_accept_index = 5
    counts_index = 7
    assigned_counts_index = 8
    '''
    
    # This will be inefficient, because I will open each file twice...
    for f in master_fnames_list:
        with open(output_path, 'wb') as writefile:
            writer = csv.writer(writefile, delimiter='\t')
            write_headers_all_samples(sample_dir_names_list, main_dir, 
                                      chromo, f, writer)

def main():
    main_dir = sys.argv[1]
    sample_dir_names_csv = sys.argv[2]
    output_path = sys.argv[3]
    
    # Create list of sample directory names.
    sample_dir_names_list = sample_dir_names_csv.split(',')
    # Create list of chromosome names corresponding to folders within sample dir
    chromosome_list = ['chr1', 'chr2']
    # Look at one chromo first for testing. REMOVE THIS COMMENT LATER.
    chromo = chromosome_list[0]
    # Find intersection of all file names in all the samples.
    master_fnames_list = get_intersection_of_files(main_dir, 
                                                   sample_dir_names_list, 
                                                   chromo)
    combine_miso_outputs_across_samples(main_dir, sample_dir_names_list,
                                        chromo, master_fnames_list,
                                        output_path)
    

if __name__ == '__main__':
    main()