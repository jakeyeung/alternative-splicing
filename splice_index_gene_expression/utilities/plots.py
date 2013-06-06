'''
Created on 2013-06-05

@author: jyeung
'''


try: 
    from igraph import Graph as graph
except ImportError:
    print('igraph not found, ignoring...')
import matplotlib.pyplot as plt
from utilities import igraph_utilities, interactive_annotations


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
    vertex_degree_tuplelist = igraph_utilities.get_top_k_hubs(network, 3)
    # Make into one large tuple to make it easier on matplotlib
    text_tuple = tuple([i for i in vertex_degree_tuplelist])
    ax.text(0.5*max(xs), max(ys), 'Top Hubs:\n%s\n%s\n%s' %text_tuple)
    plt.show()
    
def plot_subnetwork_expression(x_vector, y_vector, 
                               bubble_size, color_vector,
                               ordered_labels, 
                               bubble_annotations,
                               xlabel, 
                               ylabel,
                               output_fullpath):
    '''
    Create bubble plot
    '''
    # Get set of color_vector
    color_set = list(set(color_vector))    # Two colors only :(
    # Plot data
    ax = plt.subplot(111)
    ax.scatter(x_vector, y_vector, s=bubble_size, linewidth=2, 
               edgecolor='w', c=color_vector, alpha=0.75)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    for x, y, lab in zip(x_vector, y_vector, ordered_labels):
        ax.text(x, y, lab, size=11, horizontalalignment='right', 
                verticalalignment='top')
    # Add annotate properties
    af =  interactive_annotations.AnnoteFinder(x_vector, y_vector, 
                                               bubble_annotations)
    plt.connect('button_press_event', af)
    # Add Legend for Each Color
    p1 = plt.Circle((0, 0), radius=bubble_size[0], color=color_set[0])
    p2 = plt.Circle((0, 0), radius=bubble_size[0], color=color_set[1])
    plt.legend([p1, p2], ['No AS Detected', 'AS Detected'])
    
    
    # Draw 45 degree line on graph. 
    xmin, xmax, ymin, ymax = plt.axis()
    vector45deg = [int(min(xmin, ymin)), int(max(xmax, ymax))]    # Spans axis limits.
    ax.plot(vector45deg, vector45deg, '--')
    # Plot with maximized window
    fig = plt.gcf()
    fig.set_size_inches(18.5,10.5)
    # plt.savefig(output_fullpath)
    plt.show()
    