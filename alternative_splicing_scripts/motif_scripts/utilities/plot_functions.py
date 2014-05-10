'''
Created on 2014-01-18

@author: jyeung
'''

import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import numpy as np
import sys
import collections

def plot_hline_segments(starts, stops, ypos, colors, labels,
                        xlim=[-10, 450]):
    '''
    Plot horizonatl line segments, start and stops represent
    lists containing xmin and xmax of line segment.
    ypos is list of ypositions
    colors is list of colors
    xlim is list of length 2 consisting of [xmin, xmax]
    '''
    for start, stop, y, label, color in \
        zip(starts, stops, ypos, labels, colors):
        plt.hlines(y=y, xmin=start, xmax=stop, label=label, color=color)
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = collections.OrderedDict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    plt.xlim(xlim)
    plt.show()

def plot_histogram(values_list, n_bins, mytitle, mylabel):
    '''
    Given list of values, plot histogram.
    '''
    plt.hist(values_list, n_bins, histtype='step', 
             stacked=False, fill=True, alpha=0.5,
             label=mylabel)
    plt.title(mytitle)

def plot_density(values_lists, mytitle='mytitle', labels_lists=['lab1', 'lab2'],
                 colors_list = ['blue', 'green'],
                 xlabel='xlabel', ylabel='ylabel',
                 xmin=0, xmax=400, ymin=0, ymax=0.03,
                 smoothness=0.2,
                 drawvline=False,
                 legend_pos=2,
                 showplot=True):
    '''
    Given list of values, plot histogram.
    Takes as input a list of lists, and each 
    sublist will be a density plot.
    
    drawvline: option to draw dotted verticle line at x=drawvline
    '''
    # Set matplotlib font size globally
    font = {'family': 'sans',
            'sans-serif': 'Arial'}
    matplotlib.rc('font', **font)
    size=50    # fontsize
    
    for values_list, mylabel, color in zip(values_lists, labels_lists, colors_list):
        density = gaussian_kde(values_list)
        density.covariance_factor = lambda : smoothness
        density._compute_covariance()
        xs = np.linspace(xmin, xmax, 200)
        plt.plot(xs, density(xs), antialiased=True, label=mylabel, color=color)
        plt.fill_between(xs, density(xs), alpha=0.5, zorder=5, antialiased=True, color=color)
    plt.title(mytitle, fontsize=size)
    plt.ylabel(ylabel, fontsize=size)
    plt.xlabel(xlabel, fontsize=size)
    plt.xticks(fontsize=size)
    plt.yticks(fontsize=size)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.legend(loc=legend_pos, prop={'size': size/2.})
    # Draw dotted vertical line, optional
    if drawvline is not False:
        plt.axvline(x=drawvline, color='black', linestyle='dashed')
    if showplot:
        plt.show()
        return None
    else:
        return None
    
def plot_barplot(values_list, mytitle, mylabels, ylabel, 
                 mytext1, mytext2, mytext3, 
                 ymin, ymax, width=0.35):
    '''
    Plot barplot. Works for two bars at the moment.
    
    Even creates a nice p-value comparison between the two bars
    '''
    # Set matplotlib font size globally
    font = {'family': 'sans',
            'sans-serif': 'Arial'}
    matplotlib.rc('font', **font)
    
    ind = np.arange(len(values_list)) + width
    
    fig, ax = plt.subplots()
    ax.p1 = plt.bar(ind, values_list, width=width, facecolor='#777777', 
                    ecolor='black')
    # Fill whitespace in the margins by adjusting subplot
    fig.subplots_adjust(bottom=0.08)
    fig.subplots_adjust(left=0.08)
    fig.subplots_adjust(right=0.99)
    fig.subplots_adjust(top=0.92)
    
    # Add grid
    ax.grid(True)
    
    # Set labels and title
    size=50
    plt.xticks(ind + width/2, mylabels, fontsize=size)
    plt.xlim(0, ind[-1] + width*2)
    plt.ylabel(ylabel, fontsize=size)
    plt.title(mytitle, fontsize=size)
    plt.ylim(ymin, ymax)
    plt.yticks(fontsize=size)
    plt.text(ind[0] + width/2, values_list[0]*1.05, mytext1, 
             ha='center', va='bottom', fontsize=size)
    plt.text(ind[1] + width/2, values_list[1]*1.05, mytext2,
             ha='center', va='bottom', fontsize=size)

    # Function for creating p-value comparison
    def label_diff(i, j, text, X, Y, width, fsize=size):
        #x = (X[i]+X[j])/2
        y = 1.1*max(Y[i], Y[j])
        #dx = abs(X[i]-X[j])
        props = {'connectionstyle':'bar','arrowstyle':'-',\
                     'shrinkA':20,'shrinkB':20,'lw':6}
        ax.annotate(text, xy=(X[i] + width, y + 0.27), zorder=10, fontsize=fsize)
        ax.annotate('', xy=(X[i] + width, y), xytext=(X[j] + width, y), arrowprops=props)
    
    # Call the function
    label_diff(0, 1, mytext3, ind, values_list, width=width/2, fsize=size)
    
    