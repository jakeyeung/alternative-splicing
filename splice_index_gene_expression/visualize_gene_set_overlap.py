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
from matplotlib_venn import venn3
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
        print('Three OptDis output partial paths'\
              '(relative to output directory) must be specified in command line.')
        sys.exit()
    gene_exprs_only_partialpath = sys.argv[1]
    probe_and_gene_partialpath = sys.argv[2]
    si_and_gene_partialpath = sys.argv[3]
    
    gene_exprs_only_fullpath = os.path.join(mydirs.outputdir, 
                                            gene_exprs_only_partialpath)
    probe_and_gene_fullpath = os.path.join(mydirs.outputdir, 
                                           probe_and_gene_partialpath)
    si_and_gene_fullpath = os.path.join(mydirs.outputdir, 
                                        si_and_gene_partialpath)
    
    # Get gene list from all three optdis results. 
    gene_exprs_only_glist = get_gene_list_from_optdis_results(gene_exprs_only_fullpath, 
                                                              headers=False)
    probe_and_gene_glist = get_gene_list_from_optdis_results(probe_and_gene_fullpath, 
                                                             headers=False)
    si_and_gene_glist = get_gene_list_from_optdis_results(si_and_gene_fullpath, 
                                                          headers=False)
    
    # Find two-way intersection between the three gene lists.
    gene_only_gene_and_probe_glist = [gene for gene in gene_exprs_only_glist\
                                       if gene in probe_and_gene_glist]
    gene_only_gene_and_si_glist = [gene for gene in gene_exprs_only_glist\
                                   if gene in si_and_gene_glist]
    probe_gene_and_si_gene_glist = [gene for gene in probe_and_gene_glist\
                                    if gene in si_and_gene_glist]
    
    # Find three way intersection between three gene lists.
    gene_probe_si_glist = [gene for gene in gene_exprs_only_glist\
                           if gene in probe_gene_and_si_gene_glist]
    
    print len(gene_exprs_only_glist), len(probe_and_gene_glist), len(si_and_gene_glist)
    print len(gene_only_gene_and_probe_glist), len(gene_only_gene_and_si_glist), len(probe_gene_and_si_gene_glist)
    print len(gene_probe_si_glist)
    
    # Arrange counts in the form of (Abc, aBc, ABc, abC, AbC, aBC, ABC)
    Abc = len(gene_exprs_only_glist)
    aBc = len(probe_and_gene_glist)
    ABc = len(gene_only_gene_and_probe_glist)
    abC = len(si_and_gene_glist)
    AbC = len(gene_only_gene_and_si_glist)
    aBC = len(probe_gene_and_si_gene_glist)
    ABC = len(gene_probe_si_glist)
    
    # Plot Venn Diagram
    venn3(subsets = (Abc, aBc, ABc, abC, abC, AbC, ABC), 
          set_labels = ('Gene Exprs Only', 'Probe and Gene Exprs', 
                        'SI and Gene Exprs'))
    plt.show()
    
    
    