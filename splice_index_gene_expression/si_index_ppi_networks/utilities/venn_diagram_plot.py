'''
jakeyeung
July 8 2013
plot venn diagrams.
'''


from matplotlib_venn import venn2, venn3
from matplotlib import pyplot as plt


def plot_two_set_venn_diagram(Ab, aB, AB, set_labels):
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