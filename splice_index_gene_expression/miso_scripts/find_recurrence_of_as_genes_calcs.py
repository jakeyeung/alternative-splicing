'''
Created on 2013-09-09

@author: jyeung

After running find_recurrence_of_as_genes.py, open its output file
and do some analytics and maybe make a pretty graph.
'''


import sys
import csv

def write_meds_to_file(inclusion_list, exclusion_list, output_path):
    writecount = 0
    with open(output_path, 'wb') as writefile:
        writer = csv.writer(writefile, delimiter='\t')
        header = ['psi', 'group', 'inclusion_or_exclusion']
        writer.writerow(header)
        for tup, i_or_e in zip(inclusion_list+exclusion_list, 
                               ['inclusion']*len(inclusion_list) + \
                               ['exclusion']*len(exclusion_list)):
            g1 = tup[0]
            g2 = tup[2]
            
            '''
            print i_or_e
            print tup[1]
            print g1
            print g2
            print s_interpolated
            raw_input()
            '''
            
            if i_or_e == 'inclusion':
                g1_row = ['0', 'PC', i_or_e]
                g2_row = ['1', 'NEPC', i_or_e]
                s_interpolated = (tup[1] - g1) / (g2 - g1)
                samp_row = [str(s_interpolated), 'sample', i_or_e]
            else:
                g1_row = ['1', 'PC', i_or_e]
                g2_row = ['0', 'NEPC', i_or_e]
                s_interpolated = (tup[1] - g2) / (g1 - g2)
                samp_row = [str(s_interpolated), 'sample', i_or_e]
            '''
            g1_row = [tup[0], 'g1', i_or_e]
            g2_row = [tup[2], 'g2', i_or_e]
            samp_row = [tup[1], 's', i_or_e]
            '''
            for row in [g1_row, g2_row, samp_row]:
                writer.writerow(row)
                writecount += 1
    print('%s rows written to file: %s' %(writecount, output_path))

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
    print inclusion_list
    print exclusion_list
    
    g1_exclusion = []
    g2_inclusion = []
    g2_exclusion = []
    g1_inclusion = []
    samp_in_g1 = []
    samp_in_g2 = []
    for g1, s, g2 in zip(group1, samp, group2):
        if g2 - g1 > 0:    # inclusion in g2
            g1_exclusion.append(g1)
            g2_inclusion.append(g2)
            if g2 - s < s - g1:    # Then it is closer to g2
                samp_in_g2.append(samp)
            elif g2 - s > s - g1:    # closer to g1
                samp_in_g1.append(samp)
            else:
                print('They are equal? g1: %s, g2: %s, samp: %s' %(g1, g2, s))
                sys.exit()
        elif g2 - g1 < 0:    # exclusion in g2
            g1_inclusion.append(g1)
            g2_exclusion.append(g2)
            if g1 - s < s - g2:    # Then it is closer to g1
                samp_in_g1.append(samp)
            elif g1 - s > s - g2:    # closer to g2.
                samp_in_g2.append(samp)
    print g1_exclusion, g1_inclusion
    print g2_exclusion, g2_inclusion
    print len(samp_in_g1), len(samp_in_g2)
    
    return inclusion_list, exclusion_list
    '''
    dist_to_g1 = 0
    dist_to_g2 = 0
    for g1, g2, s in zip(group1, group2, samp):
        dist_to_g1 += abs(g1 - s)
        dist_to_g2 += abs(g2 - s)
    print dist_to_g1
    print dist_to_g2
    '''

def read_info(input_path):
    with open(input_path, 'rb') as readfile:
        reader = csv.reader(readfile, delimiter='\t')
        header = reader.next()
        samp_psi_tup_list = []
        for row in reader:
            try:
                sample_psi = float(row[header.index('sample_psi_mean')])
            except IndexError:
                continue
            event = row[header.index('event')]
            gsymbol = row[header.index('gsymbol')]
            group = row[header.index('group')].split(',')
            psi_meds = row[header.index('psi_median')].split(',')
            group1_psi_avg, group2_psi_avg = get_group_psi_med(group, psi_meds)
            samp_psi_tup_list.append((group1_psi_avg, sample_psi, 
                                      group2_psi_avg))
        inclusion_list, exclusion_list = \
            get_subtype_inclusion_or_exclusion(samp_psi_tup_list)
        return inclusion_list, exclusion_list
        

def main(input_path, out_path):
    inclusion_list, exclusion_list = read_info(input_path)
    write_meds_to_file(inclusion_list, exclusion_list, out_path)
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('input file containing recurrence of genes must be specified'\
              ' in command line.')
        sys.exit()
    input_path = sys.argv[1]
    out_path = sys.argv[2]
    main(input_path, out_path)