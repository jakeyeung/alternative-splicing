'''
Created on 2013-12-17

@author: jyeung
'''

import matplotlib.pyplot as plt
import interactive_plot_utils

def plot_bubble_plot(x_vector, y_vector, 
                     color_vector,
                     bubble_annotations,
                     xlabel, 
                     ylabel,
                     title,
                     legend,
                     saveplot=False,
                     output_fullpath=None,
                     annotated_gene_list=None):
    '''
    Create bubble plot
    '''
    
    # Plot data
    ax = plt.subplot(111)
    
    for x, y, col in zip(x_vector, y_vector, color_vector):
        ax.scatter(x, y, linewidth=2, 
           edgecolor='w', s=150, c=col, alpha=0.4)
        
    ax.set_xlabel(xlabel, size=33)
    ax.set_ylabel(ylabel, size=33)
    ax.set_ylim([-10, 10])
    ax.set_xlim([-15, 15])
    
    if annotated_gene_list is None:
        # Add annotate properties
        af = interactive_plot_utils.AnnoteFinder(x_vector, y_vector, bubble_annotations)
        plt.connect('button_press_event', af)
    else:
        # manual annotations:
        for x, y, annote in zip(x_vector, y_vector, bubble_annotations):
            if annote in annotated_gene_list:
                t = ax.text(x,y, "%s"%(annote), horizontalalignment='left', 
                            verticalalignment='bottom', size=15)
    

    plt.title(title, size=33)
    plt.tick_params(labelsize=33)
    
    # Draw x and y intercept
    xmin, xmax, ymin, ymax = plt.axis()
    # xint
    ax.plot([0, 0], [ymin, ymax], '--')
    # y-int
    ax.plot([xmin, xmax], [0, 0], '--')

    # Add legend
    # Add Legend for Each Color
    color_set = list(set(color_vector))
    if len(color_set) != 1:
        p1 = plt.Circle((0, 0), radius=2, color='b')
        p2 = plt.Circle((0, 0), radius=2, color='r')
        plt.legend([p1, p2], legend, prop={'size': 33})
    else:
        p1 = plt.Circle((0, 0), radius=2, color='r')
        plt.legend([p1], legend, prop={'size': 33})
    
    # Plot with maximized window
    fig = plt.gcf()
    fig.set_size_inches(18.5,10.5)
    if saveplot == False:
        plt.show()
    elif saveplot == True:
        plt.savefig(output_fullpath)
        print('Plot saved to %s' %output_fullpath)
        
def plot_bubble_chart(x_vector, y_vector, annotes):
    plt.scatter(x_vector, y_vector)
    af = interactive_plot_utils.AnnoteFinder(x_vector, y_vector, annotes)
    plt.connect('button_press_event', af)

    plt.show()