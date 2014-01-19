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
                 mytext1, mytext2, ymin, ymax, width=0.35):
    '''
    Plot barplot. Works for two bars at the moment.
    '''
    ind = np.arange(len(values_list)) + width
    plt.bar(ind, values_list, width=width, facecolor='#777777', 
            ecolor='black')
    plt.xticks(ind + width/2, mylabels)
    plt.xlim(0, ind[-1] + width*2)
    plt.ylabel(ylabel)
    plt.title(mytitle)
    plt.ylim(ymin, ymax)
    plt.text(ind[0] + width/2, values_list[0]*1.05, mytext1, 
             ha='center', va='bottom')
    plt.text(ind[1] + width/2, values_list[1]*1.05, mytext2,
             ha='center', va='bottom')
    