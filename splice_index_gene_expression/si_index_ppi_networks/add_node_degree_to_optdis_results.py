'''
Created on 2013-06-05

@author: jyeung

From OptDis output and PPI network, append OptDis output 
so that the genes have a node degree corresponding to it. 
'''


import sys
import os
from utilities import set_directories, igraph_utilities, read_write_gene_info

# Constants
inputdir = 'input'
outputdir = 'output'
outputheaders = ['gene_symbol', 'subnetwork', 'node_degree']
# Set directories
mydirs = set_directories.mydirs(inputdir, outputdir)


def read_gene_write_degrees(readpath, writepath, vertex_degree_dict):
    '''
    From optdis outpath as readpath and empty write text file, we will
    read write a new text file containing the optdis outpath information and
    in addition each gene's corresponding node degree. 
    '''
    readwrite_obj = read_write_gene_info.read_write_gene_info(readpath, 
                                                              writepath, 
                                                              header=False)
    with readwrite_obj:
        # Write columns
        readwrite_obj.writenext(outputheaders)
        while True:
            try:
                row = readwrite_obj.readnext()
            except StopIteration:
                print('%s rows read, %s rows written, breaking loop' \
                      %(readwrite_obj.readrowcount, readwrite_obj.writerowcount))
                break
            # Find node degree of the gene
            gene = row[0]    # First column contains gene name.
            # Append node degree to read row. 
            row.append(vertex_degree_dict[gene])
            # Write rows
            readwrite_obj.writenext(row)
    return None

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Optdis Output, Gene Network, Output filepath must be provided in command line.' \
              '\n Paths are relative to the output directory.')
        sys.exit()
    optdis_output_partialpath = sys.argv[1]
    gene_network_partialpath = sys.argv[2]
    output_partialpath = sys.argv[3]
    
    optdis_fullpath = os.path.join(mydirs.outputdir, optdis_output_partialpath)
    gene_network_fullpath = os.path.join(mydirs.outputdir, gene_network_partialpath)
    output_fullpath = os.path.join(mydirs.outputdir, output_partialpath)
    
    g = igraph_utilities.igraph_object(gene_network_fullpath, weights=False, directed=False)
    
    vertex_degree_dict = g.get_vertex_degree_dict()
    
    read_gene_write_degrees(optdis_fullpath, output_fullpath, vertex_degree_dict)
    