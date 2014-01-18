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
    
def plot_barplot(values_list, mytitle, mylabels, ylabel, width=0.35):
    '''
    Plot barplot
    '''
    # ind = np.arange(len(values_list))    # x locations for groups
    
    _, ax = plt.subplots()
    ind = [1, 10, 15, 20]
    values_list.insert(0, sys.float_info.epsilon)
    values_list.append(sys.float_info.epsilon)
    print ind, values_list
    plt.bar(ind, values_list)
    ax.set_xticklabels(mylabels)
    ax.set_ylabel(ylabel)
    ax.set_title(mytitle)
    # ax.legend(values_list, mylabels)