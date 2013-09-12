'''
Created on 2013-06-22

@author: jyeung
Plots venn diagrams of gene lists from two optdis (as modules) results.
REQUIRES matplotlib_venn library and matplotlib.
'''


import sys
from utilities import plots
from find_different_modules import get_modules_from_file, \
    find_super_unique_modules, find_unique_and_common_modules


# Set constants for plotting.
barplot_categories = ['Non-Unique Modules', 'Unique Modules']
barplot_events = ['Gene Expression Only', 'Probe and Gene Expression']
barplot_barwidths = 0.15
color_vector = ['red', 'blue']
xlabel = ''
ylabel = 'Number of Modules'
title = 'Unique Modules in OptDis Outputs'


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('Two optdis output fullpaths (for reading), number of modules, '\
              'and venn diagram labels (CSV format) must be provided '\
              'in commandline.')
        sys.exit()
    optdis_output1 = sys.argv[1]
    optdis_output2 = sys.argv[2]
    try:
        number_of_modules = int(sys.argv[3])
    except ValueError:
        sys.exit('Number of modules must be integer, not %s' %number_of_modules)
    output_fullpath = sys.argv[4]
    
    # Get modules from file.
    modules1 = get_modules_from_file(optdis_output1, number_of_modules)
    modules2 = get_modules_from_file(optdis_output2, number_of_modules)
    
    # Find unique modules for modules 1 and 2.
    unique_modules1, common_modules1 = find_unique_and_common_modules(modules2, modules1)
    unique_modules2, common_modules2 = find_unique_and_common_modules(modules1, modules2)
    
    # Are there TRULY unique modules?
    super_uniques2 = find_super_unique_modules(unique_modules2, modules1)
    super_uniques1 = find_super_unique_modules(unique_modules1, modules2)
    
    # Plot stacked bar plot.
    non_super_unique_count1 = len(modules1) - len(super_uniques1)
    non_super_unique_count2 = len(modules2) - len(super_uniques2)
    super_unique_count1 = len(super_uniques1)
    super_unique_count2 = len(super_uniques2)
    
    plots.plot_stacked_bar([non_super_unique_count1, non_super_unique_count2], 
                           [super_unique_count1, super_unique_count2], 
                           barplot_categories, barplot_events, barplot_barwidths, 
                           color_vector, xlabel, ylabel, title, output_fullpath)
    
    