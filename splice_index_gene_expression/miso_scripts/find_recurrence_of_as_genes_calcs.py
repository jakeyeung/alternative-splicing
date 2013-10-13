'''
Created on 2013-09-09

@author: jyeung

After running find_recurrence_of_as_genes.py, open its output file
and do some analytics and maybe make a pretty graph.

What kind of analytics? psi_avg distances
'''


import sys
import csv

def write_meds_to_file(inclusion_list, exclusion_list, output_path):
    '''
    Given inclusion and exclusion list, calculate
    the distance with respect to PCa (g1) and 
    write it to file.
    '''
    writecount = 0
    distance_list = []
    with open(output_path, 'wb') as writefile:
        writer = csv.writer(writefile, delimiter='\t')
        header = ['psi', 'group', 'inclusion_or_exclusion']
        writer.writerow(header)
        for tup, i_or_e in zip(inclusion_list+exclusion_list, 
                               ['inclusion']*len(inclusion_list) + \
                               ['exclusion']*len(exclusion_list)):
            g1 = tup[0]
            g2 = tup[2]
            
            if i_or_e == 'inclusion':
                g1_row = ['0', 'PC', i_or_e]
                g2_row = ['1', 'NEPC', i_or_e]
                s_interpolated = (tup[1] - g1) / (g2 - g1)
                # Calculate distance from PC or group 1
                distance_list.append(float(s_interpolated) - 0)
                samp_row = [str(s_interpolated), 'sample', i_or_e]
            else:
                g1_row = ['1', 'PC', i_or_e]
                g2_row = ['0', 'NEPC', i_or_e]
                s_interpolated = (tup[1] - g2) / (g1 - g2)
                # Calculate distance from PC or group 1
                distance_list.append(1 - float(s_interpolated))
                samp_row = [str(s_interpolated), 'sample', i_or_e]
            for row in [g1_row, g2_row, samp_row]:
                writer.writerow(row)
                writecount += 1
        # Calculate mean of distance values.
        dist_mean = float(sum(distance_list)) / len(distance_list)
    print('%s rows written to file: %s' %(writecount, output_path))
    return dist_mean

def get_group_psi_med(group, psi_med_list):
    '''
    Get psi medians for two groups.
    '''
    group1_psi_med = []
    group2_psi_med = []
    for g, psi_med in zip(group, psi_med_list):
        if g == '1':
            group1_psi_med.append(float(psi_med))
        elif g == '2':
            group2_psi_med.append(float(psi_med))
    group1_psi_avg = float(sum(group1_psi_med)) / len(group1_psi_med)
    group2_psi_avg = float(sum(group2_psi_med)) / len(group2_psi_med)
    return group1_psi_avg, group2_psi_avg

def get_subtype_inclusion_or_exclusion(samp_psi_tup_list):
    '''
    From a tuple list of form:
    (group1_psi_avg(PCa), sample_psi, group2_psi_avg(NEPC)),
    determine if the sample_psi belongs to inclusion or exclusion.
    
    We assume that:
    g1: PCa
    g2: NEPC
    
    If g2 - g1 > 0, it is INCLUSION with respect to g2 (NEPC)
    If g2 - g1 < 0, it is EXCLUSION with respect to g2 (NEPC)
    '''
    group1 = [i[0] for i in samp_psi_tup_list]
    samp = [i[1] for i in samp_psi_tup_list]
    group2 = [i[2] for i in samp_psi_tup_list]
    
    inclusion_list = []
    exclusion_list = []
    for g1, g2, s in zip(group1, group2, samp):
        if g2 - g1 > 0:    # inclusion in g2
            inclusion_list.append((g1, s, g2))
        if g2 - g1 < 0:
            exclusion_list.append((g1, s, g2))    
    return inclusion_list, exclusion_list

def read_info(input_path):
    '''
    From miso summary of a particular sample (containing DS events of
    interest (obtained from find_recurrence_of_as_genes.py.
    
    First, loops through each row, creating a tuple list of form:
    (group1_psi_avg(PCa), sample_psi, group2_psi_avg(NEPC)
    
    Then separate the list to whether it is inclusion or exclusion
    with respect to NEPC.
    
    Outputs:
    inclusion_list:
        A list containing tuples of the form:
        (group1_psi_avg(PCa), sample_psi, group2_psi_avg(NEPC)
        where group2_psi_avg - group1_psi_avg > 0
    exclusion_list:
        Tuple list of form:
        (group1_psi_avg(PCa), sample_psi, group2_psi_avg(NEPC)
        where group2_psi_avg - group1_psi-avg < 0
    '''
    with open(input_path, 'rb') as readfile:
        reader = csv.reader(readfile, delimiter='\t')
        header = reader.next()
        samp_psi_tup_list = []
        for row in reader:
            try:
                sample_psi = float(row[header.index('sample_psi_mean')])
            except IndexError:
                continue
            group = row[header.index('group')].split(',')
            psi_meds = row[header.index('psi_median')].split(',')
            group1_psi_avg, group2_psi_avg = get_group_psi_med(group, psi_meds)
            samp_psi_tup_list.append((group1_psi_avg, sample_psi, 
                                      group2_psi_avg))
        inclusion_list, exclusion_list = \
            get_subtype_inclusion_or_exclusion(samp_psi_tup_list)
        return inclusion_list, exclusion_list
        
def write_dist_mean_to_file(dist_mean, write_path):
    with open(write_path, 'wb') as writefile:
        mywriter = csv.writer(writefile, delimiter='\t')
        mywriter.writerow([dist_mean, write_path])
    print 'Dist mean written to file: %s' %write_path

def main(input_path, out_path, dist_mean_path):
    # Grab inclusion and exclusion list
    inclusion_list, exclusion_list = read_info(input_path)
    dist_mean = write_meds_to_file(inclusion_list, exclusion_list, out_path)
    write_dist_mean_to_file(dist_mean, dist_mean_path)
    print '%s mean distance for %s' %(dist_mean, input_path)
    
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('input file containing recurrence of genes must be specified'\
              ' in command line.')
        sys.exit()
    input_path = sys.argv[1]
    out_path = sys.argv[2]
    dist_mean_path = sys.argv[3]
    main(input_path, out_path, dist_mean_path)