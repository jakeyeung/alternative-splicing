'''
Created on 2013-09-11

@author: jyeung
Input a miso event, a list of samples you want to look at, 
then output a text file containing mean psis, 95 CIs and 
counts.
'''

import sys
import csv
import os

def get_samples_from_file(sample_path):
    '''
    TODO: this is already in utilities, please delete 
    and use new tool in utilities
    From a text file of one column containing sample names, 
    iterate rows and get the samples.
    '''
    samples = []
    with open(sample_path, 'rb') as readfile:
        reader = csv.reader(readfile, delimiter='\t')
        for row in reader:
            samples += row
    return samples

def find_event_get_info(miso_event, reader, header, samp):
    '''
    Itereate rows, if event matches miso_event, record
    miso_posterior_mean, ci_low/high, and counts.
    
    Use header to find index for the proper column names.
    '''
    # Def constants
    event_name_str = 'event_name'
    miso_mean_str = 'miso_posterior_mean'
    ci_low_str = 'ci_low'
    ci_high_str = 'ci_high'
    counts_str = 'counts'
    colnames = [event_name_str, miso_mean_str, ci_low_str, 
                      ci_high_str, counts_str]
    
    samp_event_info = {}
    
    for row in reader:
        if row[header.index(event_name_str)] != miso_event:
            pass
        else:    # Row matches event of interest, record info to dic.
            for colname in colnames:
                samp_event_info[colname] = row[header.index(colname)]
            print('%s found in sample %s.' %(miso_event, samp))
            break
    else:    # loop fell without finding event.
        for colname in colnames:
            samp_event_info[colname] = 'NA'
        print('%s not found in sample %s.' %(miso_event, samp))
    return samp_event_info

def make_path_from_samp(my_samp, jdir='summary', suffix='.miso_summary'):
    '''
    From my_samp, create a relative path that can be used to open
    the miso summary file.
    '''
    suffixed_samp = ''.join([my_samp, suffix])
    rel_path = os.path.join(my_samp, jdir, suffixed_samp)
    return rel_path

def get_event_info_from_sample(miso_event, miso_summary_path, samp):
    '''
    Given miso event, open miso_summary_path/samp file and 
    get psi vals and counts.
    '''
    # Create a filename from sample to open miso summary file.
    samp_rel_path = make_path_from_samp(samp, jdir='summary', 
                                        suffix='.miso_summary')
    with open(os.path.join(miso_summary_path, samp_rel_path), 'rb') as readfile:
        reader = csv.reader(readfile, delimiter='\t')
        header = reader.next()
        samp_event_info = find_event_get_info(miso_event, reader, header, samp)
    return samp_event_info

def write_subkeys(dic, write_obj):
    '''
    Assume dic has subdics, all containing the same values, so if we
    randomly picked a dic and looked at the subdic's values, we can
    use that to write subkeys as header.
    '''
    representative_subdic = dic[dic.keys()[0]]
    write_obj.writerow(['sample'] + representative_subdic.keys())
    return None

def write_samp_info(samp, samp_info_dic, write_obj):
    '''
    Given sample and sample info, write the info to file.
    '''
    samp_info = []
    for key in samp_info_dic:
        samp_info.append(samp_info_dic[key])
    write_obj.writerow([samp] + samp_info)

def write_to_dic_to_file(dic, output_path):
    '''
    dic contains a subdic. We want key to be rownames, subkeys to be
    column names. 
    '''
    with open(output_path, 'wb') as outpath:
        out_writer = csv.writer(outpath, delimiter='\t')
        write_subkeys(dic, out_writer)
        write_count = 0
        for samp, subdic in dic.iteritems():
            write_samp_info(samp, subdic, out_writer)
            write_count += 1
    print('%s samples written to file: %s' %(write_count, output_path))
    return write_count
    
def main(miso_event, sample_path, miso_summary_path, output_path):
    my_samples = get_samples_from_file(sample_path)
    event_info = {}
    for samp in my_samples:
        samp_event_info = get_event_info_from_sample(miso_event, 
                                                    miso_summary_path, 
                                                    samp)
        event_info[samp] = samp_event_info
    write_to_dic_to_file(event_info, output_path)

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('MISO event, .txt file of samples and output path'\
              ' must be specified in command line.')
        sys.exit()
    miso_event = sys.argv[1]
    miso_gene_name = sys.argv[2]
    sample_path = sys.argv[3]
    miso_summary_path = sys.argv[4]
    output_path = sys.argv[5]
    # append output_path with miso_event
    out_filename = ''.join([miso_gene_name, '_samples_info.txt'])
    output_path_appended = os.path.join(output_path, out_filename)
    main(miso_event, sample_path, miso_summary_path, output_path_appended)