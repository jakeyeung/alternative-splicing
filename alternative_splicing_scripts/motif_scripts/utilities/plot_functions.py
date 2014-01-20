'''
Created on 2014-01-18

@author: jyeung
'''

import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import numpy as np
import sys

def plot_histogram(values_list, n_bins, mytitle, mylabel):
    '''
    Given list of values, plot histogram.
    '''
    plt.hist(values_list, n_bins, histtype='step', 
             stacked=False, fill=True, alpha=0.5,
             label=mylabel)
    plt.title(mytitle)

def plot_density(values_list, mytitle, mylabel):
    '''
    Given list of values, plot histogram.
    '''
    density = gaussian_kde(values_list)
    xs = np.linspace(-7, 7, 200)
    plt.plot(xs, density(xs), label=mylabel)
    plt.title(mytitle)
    
def plot_barplot(values_list, mytitle, mylabels, ylabel, 
                 mytext1, mytext2, mytext3, 
                 ymin, ymax, width=0.35):
    '''
    Plot barplot. Works for two bars at the moment.
    
    Even creates a nice p-value comparison between the two bars
    '''
    ind = np.arange(len(values_list)) + width
    
    fig, ax = plt.subplots()
    ax.p1 = plt.bar(ind, values_list, width=width, facecolor='#777777', 
                    ecolor='black')
    # Fill whitespace in the margins by adjusting subplot
    fig.subplots_adjust(bottom=0.04)
    fig.subplots_adjust(left=0.05)
    fig.subplots_adjust(right=0.99)
    fig.subplots_adjust(top=0.95)
    
    # Add grid
    ax.grid(True)
    
    # Set labels and title
    plt.xticks(ind + width/2, mylabels, fontsize=25)
    plt.xlim(0, ind[-1] + width*2)
    plt.ylabel(ylabel, fontsize=25)
    plt.title(mytitle, fontsize=25)
    plt.ylim(ymin, ymax)
    plt.yticks(fontsize=25)
    plt.text(ind[0] + width/2, values_list[0]*1.05, mytext1, 
             ha='center', va='bottom', fontsize=25)
    plt.text(ind[1] + width/2, values_list[1]*1.05, mytext2,
             ha='center', va='bottom', fontsize=25)

    # Function for creating p-value comparison
    def label_diff(i, j, text, X, Y, width, fsize=25):
        #x = (X[i]+X[j])/2
        y = 1.1*max(Y[i], Y[j])
        #dx = abs(X[i]-X[j])
        props = {'connectionstyle':'bar','arrowstyle':'-',\
                     'shrinkA':20,'shrinkB':20,'lw':2}
        ax.annotate(text, xy=(X[i] + width, y + 0.27), zorder=10, fontsize=fsize)
        ax.annotate('', xy=(X[i] + width, y), xytext=(X[j] + width, y), arrowprops=props)
    
    # Call the function
    label_diff(0, 1, mytext3, ind, values_list, width=width/2, fsize=25)
    
    