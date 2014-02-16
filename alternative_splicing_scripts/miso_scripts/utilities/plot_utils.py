'''
Created on 2014-01-21

@author: jyeung
'''

from matplotlib_venn import venn2, venn3
import matplotlib
import matplotlib.pyplot as plt

def plot_three_set_venn(set1, set2, set3, 
                        adj_params_dic, mycolors=('r', 'g', 'b'), 
                        mylabels=None, title='Plot title'):
    '''
    Plot three circle venn diagram.
    
    adj_params_dic is of form:
    {labelid:(x_adj, yadj), ...}
    
    If adj_params_dic is None, then no adjustmenst would be made
    '''
    # Set matplotlib font size globally
    font = {'family': 'sans',
            'weight': 'bold',
            'size': 25}
    matplotlib.rc('font', **font)
    
    Abc = len(set1 - set2 - set3)
    aBc = len(set2 - set1 - set3)
    ABc = len(set1 & set2 - set3)
    abC = len(set3 - set1 - set2)
    AbC = len(set1 & set3 - set2)
    aBC = len(set2 & set3 - set1)
    ABC = len(set1 & set2 & set3)
    
    fig = plt.figure()
    # Fill whitespace in the margins by adjusting subplot
    fig.subplots_adjust(bottom=0.10)
    fig.subplots_adjust(left=0.12)
    fig.subplots_adjust(right=0.90)
    fig.subplots_adjust(top=0.90)
    ax = fig.add_subplot(111)
    
    p = venn3(subsets=(Abc, aBc, ABc, abC, AbC, aBC, ABC), 
              set_colors=mycolors, set_labels=mylabels)
    
    # Adjust textbased on adj_params_dic
    if adj_params_dic is not None:
        for labelid, adj_params in adj_params_dic.iteritems():
            label = p.get_label_by_id(labelid)
            label.set_x(label.get_position()[0] + adj_params[0])
            label.set_y(label.get_position()[1] + adj_params[1])
    
    plt.title(title)
        
    plt.show()
    
def plot_two_set_venn(set1, set2, mycolors=('r', 'g'), 
                      mylabels=None, title='Plot title'):
    '''
    Plot two circle venn diagram
    '''
    # Set matplotlib font size globally
    font = {'family': 'sans',
            'weight': 'bold',
            'size': 25}
    matplotlib.rc('font', **font)
    
    Ab = len(set1 - set2)
    aB = len(set2 - set1)
    AB = len(set1 & set2)
    
    fig = plt.figure()
    # Fill whitespace in the margins by adjusting subplot
    fig.subplots_adjust(bottom=0.10)
    fig.subplots_adjust(left=0.12)
    fig.subplots_adjust(right=0.90)
    fig.subplots_adjust(top=0.90)
    ax = fig.add_subplot(111)
    
    p = venn2(subsets=(Ab, aB, AB), set_colors=mycolors, set_labels=mylabels)
    
    plt.title(title) 
    plt.show()

def plot_two_set_venn_diagram(Ab, aB, AB, set_labels=None):
    '''
    Plots a two-set venn diagram, A and B.
    
    Inputs:
        Ab: counts in set A but not in B.
        AB: counts in both set A and B (intersection).
        aB: counts in set B but not in A.
        set_labels: a list of length 2 for set A (index 0) and set B (index 1).
    '''
    venn2(subsets=(Ab, Ab, AB), set_labels=set_labels)
    plt.show()
    
if __name__ == '__main__':
    #plot_two_set_venn_diagram(3, 2, 1, set_labels=['set1', 'set2'])
    venn2(subsets=(2, 2, 20))
    plt.show()