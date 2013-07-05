'''
Created on 2013-06-13

@author: jyeung
'''

import sys
import os
import csv
from scipy import stats
from append_optdis_probe_or_si import read_probe_si_file
from utilities import set_directories, plots

# Set constants
inputdir = 'input'
outputdir = 'output'

xlabel = 'Average $log_2$ Gene Expression: NEPC'
ylabel = 'Average $log_2$ Gene Expression: PC'
# title = 'Subnetworks and Their Average Gene Expressions'

AS_detected_string = '_AS_detected' # Label genes with alternative splice event.

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
    
def calculate_gene_avg_exprs(optdis_output_path,
                             exprs_dict, 
                             group1_indices, group2_indices, 
                             subnetwork_list, is_gene_exprs_only=False,
                             header=False, normalize_exprs=True):
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
            if is_gene_exprs_only == False:
                probe_or_gene = row[3]
            # Search genes expression in exprs_dict
            exprs = exprs_dict[genename]
            if normalize_exprs == True:
                exprs = stats.zscore(exprs)
            elif normalize_exprs == False:
                pass
            else:
                sys.exit('Normalize expression must be True or False.')
            # Split expressions into two groups
            group1_exprs = []
            group2_exprs = []
            [group1_exprs.append(exprs[i]) for i in group1_indices]
            [group2_exprs.append(exprs[i]) for i in group2_indices]
            # Take average gene expression for each
            group1_avg = float(sum(group1_exprs))/len(group1_exprs)
            group2_avg = float(sum(group2_exprs))/len(group2_exprs)
            # Append gene name as AS_detected if gene has AS event. 
            if is_gene_exprs_only == False:
                if probe_or_gene == 'all_probes':
                    pass
                elif probe_or_gene != 'all_probes':
                    # Add (AS_detected) to end of genename
                    genename += AS_detected_string
                else:
                    sys.exit('%s neither gene or probe' %probe_or_gene)
            # Append to dictionary:
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
    subnetwork_group_avg_exprs_dic = gene_group_avg_exprs_dic    # Initialize
    for subnetwork, gene_dic in subnetwork_group_avg_exprs_dic.iteritems():
        group1_exprs = [i[0] for i in gene_dic.values()]
        group2_exprs = [i[1] for i in gene_dic.values()]
        group1_subnetwork_avg = float(sum(group1_exprs))/len(group1_exprs)
        group2_subnetwork_avg = float(sum(group2_exprs))/len(group2_exprs)
        subnetwork_group_avg_exprs_dic[subnetwork] = (group1_subnetwork_avg,
                                                      group2_subnetwork_avg)
    return subnetwork_group_avg_exprs_dic


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('Optdis results, whether its gene exprs or not, expression data, '\
              'NEPC samplenames (CSV), PC samplenames (CSV), number of '\
              'subnetworks and output filename must be specified.' \
              '\nFilepath is relative to the output directory ')
        sys.exit()
    # Partial paths are relative to the output directory. 
    optdis_partialpath = sys.argv[1]
    is_gene_exprs_only = sys.argv[2]
    expression_partialpath = sys.argv[3]
    group1_samples_csv = sys.argv[4]
    group2_samples_csv = sys.argv[5]
    # group1_samples = ['X946_Urethra', 'X972_2_Penila', 'AB352', 'X963.Lpa.JP', 'X963.L.LN']
    # group2_samples = ['X1005', 'X890L', 'X945', 'X961']
    try:
        n_subnetworks = int(sys.argv[6])
    except ValueError:
        sys.exit('Number of subnetworks must be integer.')
    title = str(sys.argv[7])    # Words separated by underscore!
    title = title.replace('_', ' ')
    saveplot = sys.argv[8]
    plot_output_partialpath = sys.argv[9]
    
    # Handle csv values
    group1_samples = group1_samples_csv.split(',')
    group2_samples = group2_samples_csv.split(',')
    sample_list = group1_samples + group2_samples
    
    # Handle True False Inputs
    if is_gene_exprs_only == 'True':
        is_gene_exprs_only = True
    elif is_gene_exprs_only == 'False':
        is_gene_exprs_only = False
    else:
        sys.exit('is_gene_exprs_only must be True or False, not %s' %is_gene_exprs_only)
    if saveplot == 'True':
        saveplot = True
    elif saveplot == 'False':
        saveplot = False
    else:
        sys.exit('Saveplot arg must be neither True or False, not %s' %saveplot)
    
    # Append full paths
    optdis_fullpath = os.path.join(mydirs.outputdir, optdis_partialpath)
    exprs_fullpath = os.path.join(mydirs.outputdir, expression_partialpath)
    plot_output_fullpath = os.path.join(mydirs.outputdir, plot_output_partialpath)
    
    # Create sequence of numbers corresponding to subnetworks.
    subnetwork_list = [(i + 1) for i in range(0, int(n_subnetworks))]
    
    exprs_dict, sample_names = read_probe_si_file(exprs_fullpath, sample_list, scale_values=False, 
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
                                                        is_gene_exprs_only=is_gene_exprs_only,
                                                        header=True,
                                                        normalize_exprs=False)
    '''
    Calculate number of genes per subnetwork. Used for plotting.
    Note, I do a loop for all subnetworks in subnetwork_list, but notice 
    I do the same thing later to obtain x and y vectors. I am doing this twice
    because for some reason, when I run calculate_subnetwork_avg_exprs, my 
    gene_group_avg_exprs_dic changes to be equal to subnetwork_group_avg_exprs_dic.
    For some reason, python doesn't like that I have two dictionaries with identical 
    keynames. 
    '''
    n_genes_per_sn = []
    genenames_each_sn = []
    color_AS_events = []
    for sn in subnetwork_list:
        # Store gene names into string, for plotting. 
        gene_list = gene_group_avg_exprs_dic[sn].keys()
        genenames_each_sn.append('\n'.join(['%s:' %sn] + gene_list))
        # color subnetwork as 'lightblue' for subnetworks containing no 
        # alternative splice events. 'lightred' for subnetworks containing
        # alternative slpice events.
        for g in gene_list:
            if g.find(AS_detected_string) != -1:    # AS string found.
                color_AS_events.append('lightgreen')
                break
        else:
            color_AS_events.append('lightblue')
        # Store number of genes per sn
        n_genes_per_sn.append(len(gene_list)) 
    
    '''
    Now we modify our dictionary to contain subnetwork as keys, and average gene
    expression between group1 and group2 across genes in subnetwork as values. 
    Note, this renders gene_group_avg_exprs_dic useless because python doesnt support
    duplicate keys. 
    '''
    subnetwork_group_avg_exprs_dic = calculate_subnetwork_avg_exprs(gene_group_avg_exprs_dic)
    # Set up vectors of x and y for plotting.
    x = []
    y = []
    # x and y vectors as group1 and group2, respectively.
    for sn in subnetwork_list:
        x.append(subnetwork_group_avg_exprs_dic[sn][0])    # Group1 Avg Exprs
        y.append(subnetwork_group_avg_exprs_dic[sn][1])    # Group2 Avg Exprs
    
    # Plot results in bubble chart.
    # First, multiply subnetwork_list by a constant to make the bubbles appear larger
    bubble_size = [i*150 for i in n_genes_per_sn]
    plots.plot_subnetwork_expression(x, y, bubble_size, 
                                     color_AS_events,
                                     subnetwork_list, 
                                     genenames_each_sn, 
                                     xlabel, ylabel, title,
                                     saveplot=saveplot,
                                     output_fullpath=plot_output_fullpath)