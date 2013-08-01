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

def split_bed_by_inclusion_exclusion(bed_path, misobf_path, 
                                     inclusion_bed_path, exclusion_bed_path):
    '''
    Read bed and misobf files, read each row, determine if it is inclusion
    or exclusion, then write bed row to inclusion/exclusion accordingly.
    '''
    # Define constants
    psi_samp1_index = 1
    psi_samp2_index = 4
    misobf_event_index = 0
    bed_event_index = 3
    inclusion_str = 'inclusion'
    exclusion_str = 'exclusion'
    append_start = 10
    
    # Init read write object
    rw_obj = read_bed_misobf(bed_path, misobf_path, 
                             inclusion_bed_path, exclusion_bed_path)
    with rw_obj:
        # Get headers for bed and misobf
        bed_header = rw_obj.bedfile.readline()
        rw_obj.misobf_reader.next()    # We wont use its header.
        
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
                misobf_row = rw_obj.misobf_reader.next()
                rowcount += 1
            except StopIteration:
                print('%s rows interated.' %rowcount)
                break
            # Make sure the events are the same
            bed_event = bed_row[bed_event_index]
            misobf_event = misobf_row[misobf_event_index]
            if bed_event != misobf_event:
                print('Warning, bed event not equal to misobf event.')
                print('Bed event: %s\nMiso event: %s' \
                      %(bed_event, 
                        misobf_event))
                sys.exit()
            # Get info from bed and miso rows.
            samp1_mean = float(misobf_row[psi_samp1_index])    # PC
            samp2_mean = float(misobf_row[psi_samp2_index])    # NEPC
            
            # Find difference between NEPC and PC.
            psi_diff = float(samp2_mean - samp1_mean)
            
            # Decide whether to write to inclusion (writer1) or
            # exclusion (writer2).
            if psi_diff > 0:
                '''
                We assumed NEPC is samp2, PC is samp1.
                Therefore, psi_diff > 0 shows inclusion
                in the NEPC.
                '''
                rw_obj.writer1.writerow(bed_row)
            elif psi_diff < 0:
                '''
                Writer 2 is for exclusion isoforms...
                '''
                rw_obj.writer2.writerow(bed_row)
            else:
                print('Expected psi_diff to be greater/less than 0')
                print('PSI diff (samp2 - samp1): %s' %psi_diff)
                sys.exit()

def main():
    if len(sys.argv) < 3:
        print('bed_path (no exon/intron suffix), misobf_path must be '\
              'specified in command line.')
    bed_path = sys.argv[1]
    misobf_path = sys.argv[2]
    suffix_list_csv = sys.argv[3]    # _exon1,exon2,_exon3,_intron1,_intron2
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
                                             exclusion_bed_path)
            print('Inclusion file written to %s' %inclusion_bed_path)
            print('Exclusion file written to %s' %exclusion_bed_path)

if __name__ == '__main__':
    main()