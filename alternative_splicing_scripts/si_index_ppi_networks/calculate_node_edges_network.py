'''
Created on 2013-06-05

@author: jyeung
Calculates node degrees from a list of protein interactions.
Requires igraph.
'''


import sys
import os
from igraph import Graph as graph
from utilities import set_directories, plots

# Set constants
inputdir = 'input'
outputdir = 'output'

# Set directories
mydirs = set_directories.mydirs(inputdir, outputdir)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Filepath containing list of gene names must be provided in command line.' \
              '\n Filepath is relative to the output directory ')
        sys.exit()
    file_partialpath = sys.argv[1]
    file_fullpath = os.path.join(mydirs.outputdir, file_partialpath)
    
    network = graph.Read_Ncol(file_fullpath, names=True, weights=False, directed=False)
    network.delete_vertices('NA')
    
    plots.plot_degree_distribution(network)
    
    
    
