'''
Created on 2013-06-05

@author: jyeung
Read OptDis output, correlate this output with 
'''

import sys
import os
import csv
from append_optdis_probe_or_si import read_probe_si_file
from utilities import set_directories, plots

# Set constants
inputdir = 'input'
outputdir = 'output'
group1_samples = ['X946_Urethra', 'X972.2.Penila', 'AB352', 'X963.Lpa.JP', 'X963.L.LN']
group2_samples = ['X1005', 'X890L', 'X945', 'X961']

# Set directories
mydirs = set_directories.mydirs(inputdir, outputdir)

def get_group1_group2_indices(group1_samples, group2_samples, sample_names):
    '''
    From headers containing list of samples, separate the samples into two groups.
    One in group1, second in group2. These groups should correspond to the groups 
    used to run OptDis. 
    '''
    group1_indices = []
    group2_indices = []
    for s in sample_names:
        if s in group1_samples:
            group1_indices.append(sample_names.index(s))
        elif s in group2_samples:
            group2_indices.append(sample_names.index(s))
        else:
            sys.exit('%s not in group1 or group2' %s)
    return group1_indices, group2_indices
    
def calculate_gene_avg_exprs(optdis_output_path, exprs_dict, 
                                   group1_indices, group2_indices, 
                                   subnetwork_list,
                                   header=False):
    '''
    Calculate avg expression for EACH gene
    Return a dictionary containing subnetwork (key) and in values contains 
    another dict with gene as key, group1_avg_exprs and group2_avg_exprs as 
    a tuple. 
    Inputs:
        optdis_output_path: path to textfile containing a list of genes
        exprs_dict: dict for genename:exprs value pairs. 
        group1_indices: list of indices corresponding to 
            samples in group1 in exprs_dict.
        group2_indices: list of indices corresponding to
            samples in group2 in exprs_dict.
    Outputs:
    gene_group_avg_exprs_dic: subnetworks as keys, a dictionary in value.
        The nested dictionary contains genename as key, group1 and group2 avg
        expression tuple as value. 
    '''
    # Initialize output dic
    gene_group_avg_exprs_dic = {}
    for i in subnetwork_list:
        gene_group_avg_exprs_dic[i] = {}
    # Open file, read each row, store group1 and group2 avg exprs for each gene
    with open(optdis_output_path, 'rb') as optdis_file:
        optdis_reader = csv.reader(optdis_file, delimiter='\t')
        if header == True:
            optdis_reader.next()
        '''
        What this nested loop does.
        For loops through each row in optdis_reader.
        Keeps track of which subnetwork is happening
        '''
        for row in optdis_reader:
            genename = row[0]    # Gene name
            subnetwork = int(row[1])    # Subnetwork Number
            # Search genes expression in exprs_dict
            exprs = exprs_dict[genename]
            # Split expressions into two groups
            group1_exprs = []
            group2_exprs = []
            [group1_exprs.append(exprs[i]) for i in group1_indices]
            [group2_exprs.append(exprs[i]) for i in group2_indices]
            # Take average gene expression for each
            group1_avg = float(sum(group1_exprs))/len(group1_exprs)
            group2_avg = float(sum(group2_exprs))/len(group2_exprs)
            if subnetwork in gene_group_avg_exprs_dic:
                gene_group_avg_exprs_dic[subnetwork][genename] = (group1_avg, 
                                                                  group2_avg)
            else:
                sys.exit('%s not found as key in dictionary' %subnetwork)
    return gene_group_avg_exprs_dic

def calculate_subnetwork_avg_exprs(gene_group_avg_exprs_dic):
    '''
    This piece of code is run after calculate_gene_avg_exprs.
    For each subnetwork, we will calculate the avg gene expression.
    Inputs:
        gene_group_avg_exprs_dic: a nested dictionary with the first key:value 
        containing subnetwork number:gene_info. gene_info is a dictionary 
        with genename as key and a tuple of group1_avg_gene_exprs and 
        group2_avg_gene_exprs as value. 
    Output:
        Dictionary containing subnetworks as keys, tuple of 
        group1_avg_subnetwork_exprs, group2_avg_subnetwork_exprs
        as tuple.  
    '''
    for subnetwork, gene_dic in gene_group_avg_exprs_dic.iteritems():
        group1_exprs = [i[0] for i in gene_dic.values()]
        group2_exprs = [i[1] for i in gene_dic.values()]
        group1_subnetwork_avg = float(sum(group1_exprs))/len(group1_exprs)
        group2_subnetwork_avg = float(sum(group2_exprs))/len(group2_exprs)
        gene_group_avg_exprs_dic[subnetwork] = (group1_subnetwork_avg,
                                                group2_subnetwork_avg)
    return gene_group_avg_exprs_dic


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('Optdis results, expression data, and number of subnetworks must be specified.' \
              '\nFilepath is relative to the output directory ')
        sys.exit()
    # Partial paths are relative to the output directory. 
    optdis_partialpath = sys.argv[1]
    expression_partialpath = sys.argv[2]
    try:
        n_subnetworks = int(sys.argv[3])
    except ValueError:
        sys.exit('Number of subnetworks must be integer.')
    
    optdis_fullpath = os.path.join(mydirs.outputdir, optdis_partialpath)
    exprs_fullpath = os.path.join(mydirs.outputdir, expression_partialpath)
    # Create sequence of numbers corresponding to subnetworks.
    subnetwork_list = [(i + 1) for i in range(0, int(n_subnetworks))]
    
    exprs_dict, sample_names = read_probe_si_file(exprs_fullpath, scale_values=False, 
                                             gene_exprs_table=True)
    '''
    read_probe_si_file function returns dictionary with gene names as keys 
    and values as a list containing gene name as first element, expression as 
    the other elements. 
    Therefore, we need to do two things:
    1) we need to remove the first element from the list for each 
        key:value pair.
    2) Remove first element in header so we get nothing but sample names.
    '''
    for gene, l in exprs_dict.iteritems():
        exprs_dict[gene] = l[1:]
    sample_names = sample_names[1:]
    
    # Get index numbers for NEPC group and PC group.
    nepc_indices, pc_indices = get_group1_group2_indices(group1_samples, 
                                                           group2_samples, 
                                                           sample_names)
    # Get avg expression per group per gene.
    gene_group_avg_exprs_dic = calculate_gene_avg_exprs(optdis_fullpath, 
                                                        exprs_dict, 
                                                        nepc_indices, 
                                                        pc_indices,
                                                        subnetwork_list,
                                                        header=False)
    subnetwork_group_avg_exprs_dic = calculate_subnetwork_avg_exprs(gene_group_avg_exprs_dic)
    
    # Set up vectors of x and y for plotting.
    x = []
    y = []
    bubble_size = []
    
    
    
    
    
    
    
    
    
    
    
    
    