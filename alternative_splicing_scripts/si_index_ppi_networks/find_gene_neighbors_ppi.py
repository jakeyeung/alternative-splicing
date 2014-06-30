'''
Created on 2013-07-08

@author: jyeung

Use igraph module to find neighbors in a ppi network.
'''


import sys
import igraph


if __name__ == '__main__':
    f1 = sys.argv[1]
    f2 = sys.argv[2]
    g1 = igraph.Graph.Read(f1, format='ncol', names=True, weights=True, 
                          directed=False)
    g2 = igraph.Graph.Read(f2, format='ncol', names=True, weights=True,
                           directed=False)
    v_indices = g1.neighborhood('CD99')
    print v_indices
    # Check weights of v_indices
    for i in range(1, len(v_indices)):
        print g1.get_eid(v_indices[0], v_indices[i]), (v_indices[0], v_indices[i])
    
    print g1.vs[v_indices[1]]['name']
    print g1.vs[g1.es[g1.get_eid(0, 1)].source]['name']
    print g1.es[0].target
    '''
    print v_indices
    # for v_index in v_indices:
    #     print g1.vs[v_index]['name']
    # print g1.es['weight']
    print g1.es[0].attributes()
    print g1.get_edgelist()[430:440]
    print g1.vs[0]
    print g1.vs[1]
    print g1.vs[2]
    print g1.es[1]
    '''