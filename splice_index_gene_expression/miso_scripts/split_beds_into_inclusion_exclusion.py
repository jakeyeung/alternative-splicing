'''
Created on 2013-07-31

@author: jyeung

After running extract_coordinates_from_miso_bf.py, we realized 
that we should split 
'''


import sys
import csv
import re


class read_bed_misobf(object):
    '''
    Class for handling the reading of bed and misobf file as well as
    writing the output into another file. 
    '''
    def __init__(self, bed_file, misobf_file, write_fullpath1, write_fullpath2):
        '''
        read_fullpath_list: list containing files to be read.
        write_fullpath: filepath to write.
        
        Example:
        fullpath1 would be for inclusion in NEPC.
        fullpath2 would be for exclusion in NEPC.
        '''
        self.bedpath = bed_file
        self.misobfpath = misobf_file
        self.writepath1 = write_fullpath1
        self.writepath2 = write_fullpath2
        if write_fullpath1 in [bed_file, misobf_file] or \
            write_fullpath2 in [bed_file, misobf_file]:
            print('Read and write path the same, exiting for your protection...')
            sys.exit()
    
    def __enter__(self):
        '''
        Open read and write paths as read and write objects.
        '''
        self.bedfile = open(self.bedpath, 'rb')
        self.misobffile = open(self.misobfpath, 'rb')
        self.writefile1 = open(self.writepath1, 'wb')
        self.writefile2 = open(self.writepath2, 'wb')
        
        self.bed_reader = csv.reader(self.bedfile, delimiter='\t')
        self.misobf_reader = csv.reader(self.misobffile, delimiter='\t')
        self.writer1 = csv.writer(self.writefile1, delimiter='\t')
        self.writer2 = csv.writer(self.writefile2, delimiter='\t')
    
    def __exit__(self, exittype, exitvalue, exittraceback):
        self.bedfile.close()
        self.misobffile.close()
        self.writefile1.close()
        self.writefile2.close()
        
def create_bed_paths(bed_path, suffix_list):
    '''
    From a bed path without exon1/2/3 intron1/2 suffixes.
    Add them manually.
    e.g.
    bed_paths_list:
    pc_vs_nepc.bed -> pc_vs_nepc_exon1.bed
    bed_paths_inclusion_list:
    pc_vs_nepc.bed -> pc_vs_nepc_exon1_inclusion.bed
    bed_paths_exclusion_list:
    pc_vs_nepc.bed -> pc_vs_nepc_exon1_exclusion.bed
    '''
    # init paths lists.
    bed_paths_list = []
    bed_paths_inclusion_list = []
    bed_paths_exclusion_list = []
    
    # Split and append stuff. 
    bed_path_split = bed_path.split('.bed')[0]    # Take strings before .bed
    for suffix in suffix_list:
        bed_path_appended = ''.join([bed_path_split, suffix, '.bed'])
        bed_path_inclusion_appended = ''.join([bed_path_split, suffix, 
                                               '_inclusion.bed'])
        bed_path_exclusion_appended = ''.join([bed_path_split, suffix,
                                               '_exclusion.bed'])
        bed_paths_list.append(bed_path_appended)
        bed_paths_inclusion_list.append(bed_path_inclusion_appended)
        bed_paths_exclusion_list.append(bed_path_exclusion_appended)
    return bed_paths_list, bed_paths_inclusion_list, bed_paths_exclusion_list

def modify_bed_header(bed_header, str_to_append, append_start=10):
    '''
    Given bed header, convert
    track name = str1 -> track name = str1_str_to_append
    
    append_start = 10 because bed_header has row format:
        track name=str
        append_start corresponds to the string after =
    '''
    # Regex search track name =:
    # alphanumerics after first equal sign
    m = re.search(r'(=)\w+', bed_header)
    exprs = m.group()
    
    # Append string
    appended_exprs = '_'.join([exprs, str_to_append])
    
    # Define append finish
    append_end = append_start + len(m.group())
    
    bed_header_list = list(bed_header)
    bed_header_list[append_start:append_end] = list(appended_exprs)
    return ''.join(bed_header_list)

def check_bed_miso_events_match(bed_event, miso_event):
    if bed_event != miso_event:
        print('Warning, bed event not equal to miso event.')
        print('Bed event: %s\nMiso event: %s' \
              %(bed_event, 
                miso_event))
        sys.exit()
    return None

def decide_inclusion_or_exclusion_misobf(bed_row, summary_row):
    '''
    Input: summary_row, a row from either miso_bf or miso_t_test textfile.
        bed_row, a row from bedfile, to check if it is same event. 
        
    Output: strings 'inclusion' or 'exclusion'.
    '''
    output_str = ''
    psi_samp1_index = 1
    psi_samp2_index = 4
    miso_event_index = 0
    bed_event_index = 3
    inclusion_str = 'inclusion'
    exclusion_str = 'exclusion'
    
    bed_event = bed_row[bed_event_index]
    miso_event = summary_row[miso_event_index]
    check_bed_miso_events_match(bed_event, miso_event)
    
    # Get info from miso rows.
    samp1_mean = float(summary_row[psi_samp1_index])    # PC
    samp2_mean = float(summary_row[psi_samp2_index])    # NEPC
    # Find difference between NEPC and PC.
    psi_diff = float(samp2_mean - samp1_mean)
    if psi_diff > 0:
        '''
        We assumed NEPC is samp2, PC is samp1.
        Therefore, psi_diff > 0 shows inclusion
        in the NEPC.
        '''
        output_str = inclusion_str
    elif psi_diff < 0:
        '''
        psi_diff <= 0... 
        This therefore shows exclusion in NEPC
        '''
        output_str = exclusion_str
    else:
        print('psi_diff, %s, neither > or < 0.' %psi_diff)
    return output_str

def get_psi_diff_from_psi_list(group_list, psi_medians_list):
    '''
    group_list contains either '1' or '2' indicating group 1 
    or group 2. Normally group 1 is PCa and group 2 is NEPC.
    Find difference between the mean of the two groups's medians.
    
    psi_meds_list is ordered the same as group_list.
    
    psi_diff > 0 means inclusion in NEPC. 
    psi_diff < 0 means exclusion in NEPC.
    '''
    g1 = []
    g2 = []
    for g, psi_med in zip(group_list, psi_medians_list):
        if g == '1':
            g1.append(psi_med)
        elif g == '2':
            g2.append(psi_med)
        else:
            print('Expected %s to be "1" or "2"' %g)
            sys.exit()
    mean_g1 = float(sum(g1)) / len(g1)
    mean_g2 = float(sum(g2)) / len(g2)
    psi_diff = mean_g2 - mean_g1
    return psi_diff

def decide_inclusion_or_exclusion_t_test(bed_row, summary_row, header):
    '''
    Input: summary_row, a row from either miso_bf or miso_t_test textfile.
        bed_row, a row from bedfile, to check if it is same event. 
        
    Output: strings 'inclusion' or 'exclusion'.
    '''
    output_str = ''    # init
    # Define colname strings
    psi_median_str = 'psi_median'
    psi_median_index = header.index(psi_median_str)
    event_str = 'event'
    summary_event_index = header.index(event_str)
    group_str = 'group'
    group_index = header.index(group_str)
    # psi_median_index = 11
    # summary_event_index = 0
    # group_index = 4
    bed_event_index = 3    # Fourth column.
    
    # Get events for miso and bed
    bed_event = bed_row[bed_event_index]
    miso_event = summary_row[summary_event_index]
    check_bed_miso_events_match(bed_event, miso_event)
    
    # Get info from miso output, t_tested and filtered.
    group_list = summary_row[group_index].split(',')
    # Convert psi_meds to list and a float.
    psi_medians_list = [float(i) for i in \
                            summary_row[psi_median_index].split(',')]
    
    psi_diff = get_psi_diff_from_psi_list(group_list, psi_medians_list)
    if psi_diff > 0:
        output_str = 'inclusion'    # with respect to NEPC or group 2
    elif psi_diff < 0:
        output_str = 'exclusion'    # with respect to NEPC or group 2.
    else:
        print('%s neither greater than or less than 0.' %psi_diff)
    return output_str
            
def split_bed_by_inclusion_exclusion(bed_path, misobf_path, 
                                     inclusion_bed_path, exclusion_bed_path,
                                     hyp_test_type):
    '''
    Read bed and misobf files, read each row, determine if it is inclusion
    or exclusion, then write bed row to inclusion/exclusion accordingly.
    
    hyp_test_type must be either "bf" or "ttest".
    '''
    
    # Define constants
    '''
    psi_samp1_index = 1
    psi_samp2_index = 4
    misobf_event_index = 0
    bed_event_index = 3
    '''
    append_start = 10
    inclusion_str = 'inclusion'
    exclusion_str = 'exclusion'
    
    # Init read write object
    rw_obj = read_bed_misobf(bed_path, misobf_path, 
                             inclusion_bed_path, exclusion_bed_path)
    with rw_obj:
        # Get headers for bed and misobf
        bed_header = rw_obj.bedfile.readline()
        header = rw_obj.misobf_reader.next()    # Use header to get index.
        
        # Write bed header to first line of both output files.
        for append_str, writef in zip([inclusion_str, exclusion_str],
                                      [rw_obj.writefile1, 
                                       rw_obj.writefile2]):
            # Add inclusion or exclusion to track name description.
            bed_header_modded = modify_bed_header(bed_header,
                                                  append_str,
                                                  append_start=append_start)
            writef.write(bed_header_modded)
        
        #  Iterate read file rows in parallel.
        rowcount = 0
        while True:
            try:
                bed_row = rw_obj.bed_reader.next()
                miso_row = rw_obj.misobf_reader.next()
                rowcount += 1
            except StopIteration:
                print('%s rows interated.' %rowcount)
                break
            # Make sure the events are the same
            '''
            Refactored this into a function, output inclusion/exclusion.
            
            incl_or_excl_str is inclusion/exclusion with respect to
            group2 (NEPC). Inclusion means inclusion in NEPC and not in PC.
            '''
            if hyp_test_type=='bf':
                incl_or_excl_str = \
                    decide_inclusion_or_exclusion_misobf(bed_row, miso_row, 
                                                         header)
            elif hyp_test_type=='ttest':
                incl_or_excl_str = \
                    decide_inclusion_or_exclusion_t_test(bed_row, miso_row, 
                                                         header)
            # Decide whether to write to inclusion (writer1) or
            # exclusion (writer2).
            if incl_or_excl_str == 'inclusion':
                '''
                We assumed NEPC is samp2, PC is samp1.
                Therefore, psi_diff > 0 shows inclusion
                in the NEPC.
                '''
                rw_obj.writer1.writerow(bed_row)
            elif incl_or_excl_str == 'exclusion':
                '''
                Writer 2 is for exclusion isoforms...
                '''
                rw_obj.writer2.writerow(bed_row)
            else:
                print('Expected psi_diff to be greater/less than 0')
                print('incl_o_excl_str: %s' %incl_or_excl_str)
                sys.exit()

def main():
    if len(sys.argv) < 5:
        print('bed_path (no appended/exon/intron suffix), misobf_path must be '\
              'specified in command line.')
    bed_path = sys.argv[1]
    misobf_path = sys.argv[2]
    suffix_list_csv = sys.argv[3]    # _exon1,exon2,_exon3,_intron1,_intron2
    hyp_test_type = sys.argv[4]    # "bf" or "ttest"
    # inclusion_bed_path = sys.argv[3]    # Relative to NEPC
    # exclusion_bed_path = sys.argv[4]    # Relative to NEPC
    
    # Create list of read bed paths:
    suffix_list = suffix_list_csv.split(',')
    bed_paths_list, bed_paths_incl_list, bed_paths_excl_list = \
        create_bed_paths(bed_path, suffix_list)
        
    for bed_path, inclusion_bed_path, exclusion_bed_path in \
        zip(bed_paths_list, bed_paths_incl_list, bed_paths_excl_list):
            print('Splitting bed files to inclusion and exclusion.'\
                  '\nAssuming sample1_posterior_mean is PC and '\
                  'sample2_posterior_mean is NEPC...')
            print('Bed file: %s' %bed_path)
            print('Miso BF file: %s' %misobf_path)
            print('Outputs: \n%s\n%s' %(inclusion_bed_path, 
                                        exclusion_bed_path))
            split_bed_by_inclusion_exclusion(bed_path, misobf_path, 
                                             inclusion_bed_path, 
                                             exclusion_bed_path,
                                             hyp_test_type=hyp_test_type)
            print('Inclusion file written to %s' %inclusion_bed_path)
            print('Exclusion file written to %s' %exclusion_bed_path)

if __name__ == '__main__':
    main()