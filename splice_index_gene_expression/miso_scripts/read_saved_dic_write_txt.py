'''
Created on 2013-08-26

@author: jyeung

From the big script t_test_miso_output.py, we saved a big dictionary containing all the information
resulting from the t-test. We want to open this dic (in pickle file) and then read its contents and
filter out the events that do not have sufficient reads.

Things to filter:
number of total reads aligning to any isoform (or to both isoforms) >= N
Number of inclusion reads (i.e. reads supporting the first isoform) has to be 
    greater than or equal to N.
Number of exclusion reads (i.e. reads supporting the second isoform) has to 
    be greater than or equal to N.
The sum of inclusion and exclusion reads has to be greater than or equal to N.

MISO recommends contain only events with: 
(a) at least 1 inclusion read, 
(b) 1 exclusion read, such that 
(c) the sum of inclusion and exclusion reads is at least 10

These recommendations are for between two samples. We have two groups. But we will actually 
get the SUM of inclusion/exclusion reads from all samples (in each group) in order to 
increase the power.
'''


import sys
from group_miso_utils import read_pickle, create_chromo_list
from t_test_miso_output import read_pickle_write_to_file


def main():
    if len(sys.argv) < 2:
        print('Pickle file from t_test_miso_output.py and output '\
              '.txt filename must be specified in command line.')
        sys.exit()
    pickle_path = sys.argv[1]
    writefile_path = sys.argv[2]
    
    chr_list = create_chromo_list()
    
    # Read pickle file to get fnames_dic
    # This only contains filenames, no data. 
    fnames_dic = read_pickle(pickle_path)
    
    # Read and write to file. 
    read_pickle_write_to_file(writefile_path, chr_list, fnames_dic, 
                              filter_events=True)
    print('Summary file saved in: %s' %writefile_path)
    
if __name__ == '__main__':
    main()