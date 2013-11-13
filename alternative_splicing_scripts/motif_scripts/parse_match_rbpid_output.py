'''
Created on 2013-11-12

@author: jyeung
Reads the batch output from match_rbpid_to_rbpname (stored in a txtfile) and
extracts:
    # rbps, meme_params, location(with respect to SE), shuffled or not.

Example output:
G:\jyeung\projects\alternative_splicing\git\alternative-splicing\shellscripts\match_rbpid_to_rbpname_batch.sh
'''

import sys
import csv
from optparse import OptionParser

def list_begins_with_integer(jlist):
    '''
    Check if first element in list is an integer.
    '''
    try:
        int(jlist[0])
        return True
    except ValueError:
        return False
    
def create_colnames():
    '''
    Create a list, each element represents a colname.
    '''
    n_rbps_str = 'Number_of_RBPs'
    params_str = 'MEME_and_filter_parameters'
    location_str = 'Location_around_cassette'
    shuffled_str = 'Shuffled_or_Not'
    # Create list
    colnames = [n_rbps_str, params_str, location_str, shuffled_str]
    return colnames

def extract_dir_from_path(myfilename, extract_index=0):
    '''
    Takes filename, split by '/', then extracts the ith element
    from the split'd list.
    '''
    myfilename_split = myfilename.split('/')
    return myfilename_split[extract_index]

def main():
    usage = 'usage: %prog match_rbpid_file.txt output_file.txt'
    parser = OptionParser(usage=usage)
    (_, args) = parser.parse_args()
    if len(args) < 2:
        print 'Two args must be specified in command line:\n'\
            '1) batch output from match_rbpid_rbpname.py (like nohup output)\n'\
            '2) Output file\n'
        sys.exit()
    match_rbpid_file = args[0]
    output_file = args[1]
    
    # Define my indexes for extracting match_rbpid info.
    n_rbps_i = 0
    filename_i = 5
    
    # Initialize write file
    writefile = open(output_file, 'wb')
    mywriter = csv.writer(writefile, delimiter='\t')
    # Write header
    colnames = create_colnames()
    mywriter.writerow(colnames)
    
    # Read and iterate lines in match_rbpid_file
    writecount = 0
    with open(match_rbpid_file, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter=' ')  # Split by space.
        for row in myreader:
            # Check first "column" is an integer.
            # Rows not beginning with integer are ignored, since they may
            # be run-errors.
            if list_begins_with_integer(row):
                n_rbps = row[n_rbps_i]
                myfilename = row[filename_i]
                params = extract_dir_from_path(myfilename, 
                                               extract_index=9)
                location = extract_dir_from_path(myfilename, 
                                                 extract_index=11)
                shuffled_status = extract_dir_from_path(myfilename, 
                                                        extract_index=10)
                # Write the four extracted info to file.
                mywriter.writerow([n_rbps, params, location, shuffled_status])
                writecount += 1
    print '%s rows written to file: %s' %(writecount, output_file)
    writefile.close()

if __name__ == '__main__':
    main()
