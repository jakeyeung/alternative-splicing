'''
Created on 2013-06-05

@author: jyeung
'''

import sys
try: 
    from igraph import Graph as graph
except ImportError:
    print('igraph not found, ignoring...')
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
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
                               title,
                               saveplot=False,
                               output_fullpath=None):
    '''
    Create bubble plot
    '''
    
    # font = {'weight' : 'normal','size'   : 33}
    # matplotlib.rc('font', **font)
    
    # Get set of color_vector
    color_set = list(set(color_vector))    # Two colors only :(
    
    # Begin, some code I added for the case of probe and gene exprs...
    unique_modules = [9, 19, 22, 28, 32, 33, 35, 42, 44, 45, 46, 47, 49, 50]
    marker_list = []
    for l in ordered_labels:
        if l in unique_modules and len(color_set) == 2:
            marker_list.append('s')
        else:
            marker_list.append('o')
    # End, some code I added.
    
    # Plot data
    ax = plt.subplot(111)
    
    for x, y, size, col, mark in zip(x_vector, y_vector, bubble_size, color_vector, marker_list):
        ax.scatter(x, y, s=size, linewidth=2, 
           edgecolor='w', c=col, marker=mark, alpha=0.75)
        
    # ax.scatter(x_vector, y_vector, s=bubble_size, linewidth=2, 
    #            edgecolor='w', c=color_vector, marker='o', alpha=0.75)
    ax.set_xlabel(xlabel, size=33)
    ax.set_ylabel(ylabel, size=33)
    for x, y, lab in zip(x_vector, y_vector, ordered_labels):
        ax.text(x, y, lab, size=11, horizontalalignment='right', 
                verticalalignment='top')
    # Add annotate properties
    af =  interactive_annotations.AnnoteFinder(x_vector, y_vector, 
                                               bubble_annotations)
    plt.connect('button_press_event', af)
    # Add Legend for Each Color
    p1 = plt.Circle((0, 0), radius=bubble_size[0], color=color_set[0])
    if len(color_set) == 2:
        p2 = plt.Circle((0, 0), radius=bubble_size[0], color=color_set[1])
        p3 = plt.Circle((0, 0), radius=bubble_size[0], color=color_set[0])
        p4 = plt.Circle((0, 0), radius=bubble_size[0], color=color_set[0])
        plt.legend([p1, p2, p3, p4], ['No AS Detected', 'AS Detected', 'Unique Modules', 'Non-Unique Modules'], prop={'size':33})
    elif len(color_set) == 1:
        pass
    else:
        sys.exit('Only 1 or 2 colors can be specified.')
    plt.title(title, size=33)
    plt.tick_params(labelsize=33)
    
    # Draw 45 degree line on graph. 
    xmin, xmax, ymin, ymax = plt.axis()
    vector45deg = [int(min(xmin, ymin)), int(max(xmax, ymax))]    # Spans axis limits.
    ax.plot(vector45deg, vector45deg, '--')
    # Plot with maximized window
    fig = plt.gcf()
    fig.set_size_inches(18.5,10.5)
    if saveplot==False:
        plt.show()
    elif saveplot == True:
        plt.savefig(output_fullpath)
        print('Plot saved to %s' %output_fullpath)
    
def plot_stacked_bar(count1, count2,
                     barplot_categories, barplot_events,
                     barplot_barwidths,
                     color_vector, xlabel, ylabel, title,
                     output_fullpath):
    
    font = {'weight' : 'normal','size'   : 33}
    matplotlib.rc('font', **font)
    
    '''
    count1: count horizontally
    count2: count horizontally.
    Height of bar is count1[i] + count2[i]
    '''
    ind = range(0, len(count1))
    ind = [i + 0.5 for i in ind]
    
    count_list = [count1, count2]
    iteration = [count2, [0]*len(count2)]
    for count, category, color, bottom in zip(count_list, 
                                       barplot_categories, 
                                       color_vector, iteration):
        plt.bar(ind, tuple(count), width=barplot_barwidths, 
                color=color, label=category, bottom=bottom)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.xticks([(i + barplot_barwidths/2.) for i in ind], tuple(barplot_events))
    fig = plt.gcf()
    fig.set_size_inches(18.5,10.5)
    plt.savefig(output_fullpath)
    print('Plot saved to %s' %output_fullpath)
    plt.show()
    plt.clf()
    
def plot_line_graph(x, y, xlabel, ylabel, title):
    plt.figure(1)
    plt.plot(x, y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.show()
    
def plot_histogram(x, xlabel, ylabel, title, nbins=30, normed=False, 
                   facecolor='green', alpha=0.75, add_best_fit=False):
    '''
    From a list of values, create histogram. Requires specification of
    number of xlabel, ylabel, title.
    '''
    # Initialize a plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Plot histogram
    _, bins, _ = ax.hist(x, nbins, normed=normed, facecolor=facecolor, alpha=alpha)
    '''
    # hist uses np.histogram under the hood to create 'n' and 'bins'.
    # np.histogram returns the bin edges, so there will be 50 probability
    # density values in n, 51 bin edges in bins and 50 patches.  To get
    # everything lined up, we'll compute the bin centers
    '''
    bincenters = 0.5*(bins[1:]+bins[:-1])
    # add a 'best fit' line for the normal PDF
    if add_best_fit==True:
        mu, sigma = 100, 15    # Constants to define for bestfit
        y = mlab.normpdf( bincenters, mu, sigma)
        ax.plot(bincenters, y, 'r--', linewidth=1)
    elif add_best_fit==False:
        pass
    else:
        sys.exit('add_best_fit must be True or False, %s found' %add_best_fit)
    # Add axis labels
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    # ax.set_xlim(40, 160)
    # ax.set_ylim(0, 0.03)
    ax.grid(True)
    
    plt.show()