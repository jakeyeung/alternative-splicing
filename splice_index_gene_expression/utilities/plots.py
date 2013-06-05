'''
Created on 2013-06-05

@author: jyeung
'''


from igraph import Graph as graph
import matplotlib.pyplot as plt
from utilities import graph_utilities


def plot_degree_distribution(network):
    '''
    From an igraph object (initialized probably from graph.Read_Ncol(filename),
    return a plot with the degree distributions.
    
    Write the top 3 hubs in the network on the graph. 
    '''
    # Pair vertex with node degree, sort by degree. 
    vertex_degree_tuple = zip(network.vs, network.degree())
    vertex_degree_tuple = sorted(vertex_degree_tuple, 
                                 key=lambda x: x[1], 
                                 reverse=True)
    
    xs, ys = zip(*[(left, count) for left, _, count in \
                   graph.degree_distribution(network, bin_width=1).bins()])
    ax = plt.subplot(111)
    ax.plot(xs, ys)
    ax.set_xscale('log')
    ax.set_xlabel('Node Degree')
    ax.set_ylabel('Number of Genes')
    # Add text for top hubs
    vertex_degree_tuplelist = graph_utilities.get_top_k_hubs(network, 3)
    # Make into one large tuple to make it easier on matplotlib
    text_tuple = tuple([i for i in vertex_degree_tuplelist])
    ax.text(0.5*max(xs), max(ys), 'Top Hubs:\n%s\n%s\n%s' %text_tuple)
    plt.show()
    