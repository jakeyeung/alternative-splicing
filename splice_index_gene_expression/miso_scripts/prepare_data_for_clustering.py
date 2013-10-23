'''
Created on 2013-10-23

@author: jyeung
Read filtered t-test output and 
create matrix with gene names as rownames
and samples on columns.
Fill matrix with PSI median values, if 
not available, then write NA
'''

import sys
import csv
from optparse import OptionParser
import group_miso_utils


def get_row_info(row, rowheader, colname_str, csv=False):
    '''
    Grabs information from row. Outputs a list if it's CSV, string otherwise.
    Inputs:
        row: row from a file (like from t_test_output)
        rowheader: header of the row, to match with colname_str
        colname_str: column name of the info you want to extract.
        csv = whether that element is a CSV format or not. If it is
            csv, will try to convert and return as a list.
    '''
    my_info = row[rowheader.index(colname_str)]
    if csv == True:
        my_info = my_info.split(',')
    elif csv == False:
        pass
    else:
        print 'CSV must be True or False. %s found.' %csv
        sys.exit()
    return my_info
    
def init_samp_dic(mysamps):
    '''
    initialize a dictionary of form:
            {samp: 'NA'} 
    '''
    mydic = {}
    for s in mysamps:
        mydic[s] = 'NA'
    return mydic

def fill_dic(mydic, mykeys, myvalues):
    '''
    Input:
        mydic: initialized dic with {mykey: 'NA'}
        mykeys: keynames of dic
        myvalues: values of dic
    Output:
        mydic, filled int.
    '''
    for k, v in zip(mykeys, myvalues):
        try:
            mydic[k] = v
        except KeyError:
            print 'Error, %s was not initialized as a key in dictionary.' %k
    return mydic

def match_dic_to_header(mydic, header):
    '''
    Match keys of mydic to header, returning a list of 
    dic values in the order of the header.
    '''
    ordered_values = []
    for colname in header:
        ordered_values.append(mydic[colname])
    return ordered_values

def write_row_dic_to_file(row_dic, writeheader, mywriter):
    '''
    Form psi meds (found in row_dic) and gsymbol into a list that corresponds
    to the writeheader
    '''
    ordered_values = match_dic_to_header(row_dic, writeheader)
    mywriter.writerow(ordered_values)
    return None

def main():
    # Get parse options.
    usage = 'usage: %prog [options] t_test_file samp1_file samp2_file write_output_file'
    parser = OptionParser(usage=usage)
    (_, args) = parser.parse_args()
    
    if len(args) < 4:
        print('Following must be specified in command line:\n'\
              't_test_file (filtered, preferably)\n'\
              'samp1_file (group1 samps: should match sample names in t_test_file\n'\
              'samp2_file (group2 samps: should match samp names in t_test_file\n'\
              'output_file')
        sys.exit()
        
    # Define args to variables
    t_test_file = args[0]
    samp1_file = args[1]
    samp2_file = args[2]
    output_file = args[3]
    
    # Define colnames in ttest output
    gsymbol_str = 'gsymbol'
    sample_name_str = 'sample_name'
    psi_median_str = 'psi_median'
    
    # Get sample names from textfile.
    group_1_samples = \
        group_miso_utils.get_sample_names_from_file(samp1_file)
    group_2_samples = \
        group_miso_utils.get_sample_names_from_file(samp2_file)
    all_samples = group_1_samples + group_2_samples
        
    # Init my write file
    writefile = open(output_file, 'wb')
    mywriter = csv.writer(writefile, delimiter='\t')
    # Write row of columns
    # first column reserved for genenames therfore we use ['gsymbol']
    write_header = [gsymbol_str] + all_samples   
    mywriter.writerow(write_header)
    
    # Read t_test_file
    with open(t_test_file, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        header = myreader.next()
        writecount = 0
        for row in myreader:
            '''
            For each row, grab sample names and their PSI median values.
            Many rows will not have PSI medians due to low read count.
            Replace missing PSI medians with NAs.
            
            For each iteration, initialize a dictionary of form:
            {samp: 'NA'}
            fill in psi meds into a dictionary of form:
            {samp: psi_med}
            '''
            # Initialize dictionary
            row_dic = init_samp_dic(all_samples + [gsymbol_str])
            
            # Get gsymbol, row samps (as list), psi meds (as list)
            gsymbol = get_row_info(row, header, gsymbol_str, csv=False)
            row_samps = get_row_info(row, header, sample_name_str, csv=True)
            row_psi_meds = get_row_info(row, header, psi_median_str, csv=True)

            # Fill in rowsamps, psimeds, gsymbol into our initialized row_dic
            row_dic = fill_dic(row_dic, 
                               row_samps + [gsymbol_str], 
                               row_psi_meds + [gsymbol])
            
            # Write row_dic and gsymbol to file.
            write_row_dic_to_file(row_dic, 
                                  write_header, 
                                  mywriter)
            writecount += 1
    print('%s rows written to: %s' %(writecount, output_file))
    writefile.close()
            

if __name__ == '__main__':
    main()