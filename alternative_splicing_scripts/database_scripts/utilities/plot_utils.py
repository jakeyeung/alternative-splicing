'''
Created on 2014-01-13

@author: jyeung
'''

import matplotlib.pyplot as plt
import numpy as np

def plot_bar_plot(vector1, vector2, xticks_vector, ylabel, title, width=0.35):
    '''
    Plots barplot, bars stick together by xticks.
    
    Inputs:
    vector1: vector of floats1
    vector2: vector of floats2
    xticks_vector: xticks label, should match vector1 and vector2.
    ylabel: string to label y axis
    title: string to label title
    '''
    # make sure lenghts of vec1 vec2 xticks are the same
    for _, _, _ in zip(vector1, vector2, xticks_vector):
        pass
    
    ind = np.arange(len(vector1))
    
    fig, ax = plt.subplots()
    rects1 = ax.bar(ind, vector1, width, color='r')
    rects2 = ax.bar(ind+width, vector2, width, color='b')

    # add some
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xticks(ind+width)
    ax.set_xticklabels(xticks_vector)
    
    def autolabel(rects):
        # attach some text labels
        for rect in rects:
            height = rect.get_height()
            ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
                    ha='center', va='bottom')
    
    autolabel(rects1)
    autolabel(rects2)
    
    plt.xticks(rotation=70)
    plt.show()
    

def plot_uniprot_annots(incl_excl_dic):
    '''


    '''
    # Get full list of features, into one list.
    features_list = []
    incl_excl_list = []
    for incl_or_excl in incl_excl_dic:
        incl_excl_list.append(incl_or_excl)
        features_list += incl_excl_dic[incl_or_excl]
    # Remove redundancies.
    features_list = list(set(features_list))
    
    # Build bar heights
    n_features1 = []
    n_features2 = []
    for incl_excl, features_vector in zip(incl_excl_list, 
                                          [n_features1, 
                                           n_features2]):
        for feature in features_list:
            incl_excl_dic[incl_excl][feature].append(features_vector)
    
    # BEGIN PLOTTING
    # get indices
    ind = np.arange(len(features_list))
    width = 0.35
    
    fig, ax = plt.subplots()
    rects1 = ax.bar(ind, n_features1, width, color='b')
    rects2 = ax.bar(ind, n_features2, width, color='r')
    
    plt.show()
    
    