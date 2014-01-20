'''
Created on 2014-01-13

@author: jyeung
'''

import sys
import matplotlib.pyplot as plt
import numpy as np

def plot_bar_plot(vector1, vector2, xticks_vector, ylabel, title,
                  label1, label2, width=0.35, text_pos=0.2, stacked=True):
    '''
    Plots barplot, bars stick together by xticks.
    
    Inputs:
    vector1: vector of floats1
    vector2: vector of floats2
    xticks_vector: xticks label, should match vector1 and vector2.
    ylabel: string to label y axis
    title: string to label title
    
    width: width of bars
    text_pos: how high from top of rectangle shoudl text be
    '''
    # make sure lenghts of vec1 vec2 xticks are the same
    for _, _, _ in zip(vector1, vector2, xticks_vector):
        pass
    
    ind = np.arange(len(vector1))
    if stacked:
        ind2 = ind
    elif not stacked:
        ind2 = ind + width    # bars separated, not stacked
    else:
        print 'stacked must be True or False. %s found.' %stacked
        sys.exit()
    
    fig = plt.figure()
    # Fill whitespace in the margins by adjusting subplot
    fig.subplots_adjust(bottom=0.17)
    fig.subplots_adjust(left=0.04)
    fig.subplots_adjust(right=0.99)
    fig.subplots_adjust(top=0.95)
    ax = fig.add_subplot(111)
    #_, ax = plt.subplots()
    
    rects1 = ax.bar(ind, vector1, width, color='y')
    rects2 = ax.bar(ind2, vector2, width, color='r', bottom=vector1)

    # add some
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xticks(ind+width/5)
    ax.set_xticklabels(xticks_vector)
    
    # Adjust yaxis and title to 25 font
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(25)
    # Adjust xtick labels to 15 font
    for i in ax.get_xticklabels():
        i.set_fontsize(15)
    
    # add legend
    ax.legend((rects1[0], rects2[0]), (label1, label2), loc='upper left')
    
    # add grid
    ax.grid(True)
    
    # rotate ticks
    plt.xticks(rotation=70)
    plt.show()
    
    # optional: add text
    '''
    def autolabel(rects):
        # attach some text labels
        for rect in rects:
            height = rect.get_height()
            ax.text(rect.get_x()+rect.get_width()/2, text_pos + height, '%.0f'%float(height),
                    ha='center', va='bottom')
    autolabel(rects1)
    autolabel(rects2)
    '''
    return None

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
    
    