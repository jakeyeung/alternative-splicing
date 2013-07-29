'''
Created on 2013-07-29

@author: jyeung
'''


def append_multiple_bed_files(bed_files_list, output_file):
    '''
    List of bed files full paths, create an append bed file.
    '''
    with open(output_file, 'wb') as writefile:
        for bed_file in bed_files_list:
            with open(bed_file, 'rb') as readfile:
                for line in readfile:
                    writefile.write(line)
    return None