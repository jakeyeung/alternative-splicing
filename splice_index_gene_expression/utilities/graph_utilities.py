'''
Created on 2013-06-05

@author: jyeung
'''

import sys

def pair_vertex_with_degree(network, sortlist=True):
    '''
    From a graph object, return a list of tuples
    containing vertex names and its corresponding
    node degree, sorted or not. 
    '''
    # Pair vertex with node degree, sort by degree. 
    vertex_degree_tuple = zip(network.vs['name'], network.degree())
    if sortlist == True:
        vertex_degree_tuple = sorted(vertex_degree_tuple, 
                                     key=lambda x: x[1], 
                                     reverse=True)
    elif sortlist == False:
        pass
    else:
        sys.exit('sortlist must be either True or False.') 
    return vertex_degree_tuple

def get_top_k_hubs(network, k):
    '''
    From a graph object, return a list of tuples
    containing vertex names and its corresponding
    node degree for the top k genes ranked by node degree.
    '''
    # Pair vertex with node degree, sort by degree. 
    vertex_degree_tuple = zip(network.vs['name'], network.degree())
    vertex_degree_tuple = sorted(vertex_degree_tuple, 
                             key=lambda x: x[1], 
                             reverse=True)
    top_hubs = vertex_degree_tuple[:k]
    return top_hubs