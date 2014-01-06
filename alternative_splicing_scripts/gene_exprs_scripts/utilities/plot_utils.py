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
                     saveplot=False,
                     output_fullpath=None):
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

    # Add annotate properties
    af = interactive_plot_utils.AnnoteFinder(x_vector, y_vector, bubble_annotations)
    plt.connect('button_press_event', af)

    plt.title(title, size=33)
    plt.tick_params(labelsize=33)
    
    # Draw x and y intercept
    xmin, xmax, ymin, ymax = plt.axis()
    # xint
    ax.plot([0, 0], [ymin, ymax], '--')
    # y-int
    ax.plot([xmin, xmax], [0, 0], '--')
    
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