'''
Created on 2013-06-07

@author: jyeung

Plots venn diagrams of gene lists from three optids results.
REQUIRES matplotlib_venn library and matplotlib.
'''


import sys
import os
import csv
from utilities import set_directories
from matplotlib_venn import venn2, venn3
from matplotlib import pyplot as plt

# Set constants
inputdir = 'input'
outputdir = 'output'

# Set directories
mydirs = set_directories.mydirs(inputdir, outputdir)

def get_gene_list_from_optdis_results(optdis_fullpath, headers=False):
    gene_list = []
    with open(optdis_fullpath, 'rb') as optdis_file:
        optdis_reader = csv.reader(optdis_file, delimiter='\t')
        for row in optdis_reader:
            gene_list.append(row[0])    # First column is gene.
    return gene_list

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('Two or Three OptDis output partial paths'\
              '(relative to output directory) must be specified in command line.')
        sys.exit()
    total_data_sets = int(sys.argv[1])
    if total_data_sets not in [2, 3]:
        print('Total data sets must be either 2 or 3, %s found' %total_data_sets)
        sys.exit()
    optdis_output_1_partialpath = sys.argv[2]
    optdis_output_2_partialpath = sys.argv[3]
    venn_diagram_labels_csv = sys.argv[4]    # CSV values.
    if total_data_sets == 3:
        optdis_output_3_partialpath = sys.argv[5]
        optdis_output_partialpath_list = [optdis_output_1_partialpath, 
                                               optdis_output_2_partialpath, 
                                               optdis_output_3_partialpath]
    elif total_data_sets == 2:
        optdis_output_partialpath_list = [optdis_output_1_partialpath,
                                          optdis_output_2_partialpath]
    else:
        sys.exit('total_data_sets must be either 2 or 3, %s inputed.' \
                 %total_data_sets)
    
    # Parse CSV into a list for diagram labels.
    venn_diagram_labels_list = venn_diagram_labels_csv.split(',')
    
    # Get full paths from partial paths.
    optdis_output_fullpath_list = []
    for partialpath_output in optdis_output_partialpath_list:
        optdis_output_fullpath_list.append(os.path.join(mydirs.outputdir,
                                                        partialpath_output))
    # Get gene list from optdis results.
    output_genelists = []
    for fullpath_output in optdis_output_fullpath_list:
        output_genelists.append(get_gene_list_from_optdis_results(fullpath_output,
                                                                  headers=False))
    # Find two-way intersection between gene lists.
    common_genes = []
        
    if total_data_sets == 2:
        common_genes.append([gene for gene in output_genelists[0]\
                             if gene in output_genelists[1]])
    elif total_data_sets == 3:
        '''
        Itereate from -1 to len(genelist)-1, for finding the 
        common genes between -1 and 0, 0 and 1, 1 and 2. 
        Therefore, output list has index 0, 1, 2 with:
            0: common genes between first and third set.
            1: common genes between first and second set.
            2: common genes between second and third set.
        '''
        for i in [(j-1) for j in range(0, len(output_genelists))]:
            common_genes.append([gene for gene in output_genelists[i]\
                                 if gene in output_genelists[(i+1)]])
        '''
        With total_data_sets == 3, we have also the common genes between
        all three sets, which we will give to output list at index 3.
        '''
        common_genes.append([gene for gene in output_genelists[0]\
                           if gene in output_genelists[1] and output_genelists[2]])
    
    # Set up inputs for venn2 and venn3, which accepts count data. 
    if total_data_sets == 2:
        Ab = len(output_genelists[0])
        aB = len(output_genelists[1])
        AB = len(common_genes[0])
        venn2(subsets=(Ab, aB, AB), set_labels=venn_diagram_labels_list)
    elif total_data_sets == 3:
        Abc = len(output_genelists[0])
        aBc = len(output_genelists[1])
        abC = len(output_genelists[2])
        AbC = len(common_genes[0])
        ABc = len(common_genes[1])
        aBC = len(common_genes[2])
        ABC = len(common_genes[3])
        venn2(subsets = (Abc, aBc, ABc, abC, abC, AbC, ABC), 
              set_labels = venn_diagram_labels_list)
    plt.show()
    
    
        
        
    
    