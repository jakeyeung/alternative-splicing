'''
Created on 2014-01-18

@author: jyeung
'''

import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import numpy as np

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
    '''
    plt.hist(values_list, n_bins, histtype='step', 
             stacked=False, fill=True, alpha=0.5,
             label=mylabel)
    '''
    plt.title(mytitle)