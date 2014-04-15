'''
Created on 2014-04-01

@author: jyeung

Plot meme motif positions in a pretty plot

Reads as input the pkl file from meme_gerp_summary_pkl (includes GERP scores)
for each miso event, we have all matching regions containing start and end
sites. We will plot each miso event with its motif.
'''

import sys
from optparse import OptionParser
import pickle
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from motif_scripts.utilities import plot_functions, gerp_utilities
from plot_meme_motif_null_comparison import get_dic_from_pklpath

def add_rectangles(rect_start, height=0.002, length=10, color='orange'):
    '''
    Add rectangles to figure. Simulates
    a cassette exon cartoon.
    Rect start is top right corner x position of rectangle.
    Assumes top is at y=0.
    '''
    rect_start
    # add rectangles
    verts = [
    (rect_start-length, -height), # left, bottom
    (rect_start-length, 0.), # left, top
    (rect_start, 0.), # right, top
    (rect_start, -height), # right, bottom
    (0., 0.), # ignored
    ]
    codes = [Path.MOVETO,
             Path.LINETO,
             Path.LINETO,
             Path.LINETO,
             Path.CLOSEPOLY,
             ]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor=color, lw=1)
    return patch

def main():
    usage = 'usage: %prog meme_gerp_summary_pkl\n'\
        'Requires one argument:\n'\
        '1) pkl file from summarize_meme_results'
    parser = OptionParser(usage=usage)
    parser.add_option('-r', '--plot_raw_locations', dest='raw_locations',
                      default=False,
                      help='Boolean value. True=horizontal '\
                        'line segment plot. False=density plot') 
    (options, args) = parser.parse_args()
    
    if len(args) != 1:
        print 'Requires 1 argument to be specified in commandline'
        print usage
        sys.exit()
    pklpath = args[0]
    
    if options.raw_locations in ['True', 'TRUE', True]:
        raw_locations = True
        print 'Plotting raw locations...'
    elif options.raw_locations in ['False', 'FALSE', False]:
        raw_locations = False
        print 'Plotting density plot...'
    else:
        print '--plot_raw_locations option must be '\
            'True or False. %s found.' %options.raw_locations
        sys.exit()
    
    # get dics from pkl 
    meme_dic = get_dic_from_pklpath(pklpath)
    print meme_dic
    
    event_count = 0    # used as y-axis locater...
    # init offsetters
    offset_length = 100
    offsets = [0, 110, 230, 340, 450]
    '''
    Set Motif 1 to #CC6666, Motif 2 to #33CCCC Motif 3 to "green"
    '''
    
    plot_settings_dic = {'intron_1_5p': {'offset': offsets[0], 
                                         'color': ['#CC6666', '#CC6666', 'cyan']},
                         'intron_1_3p': {'offset': offsets[1], 
                                         'color': ['green']},
                         'intron_2_5p': {'offset': offsets[2], 
                                         'color': ['red']},
                         'intron_2_3p': {'offset': offsets[3], 
                                         'color': ['#33CCCC', 'black']}}
                        
    # collect plot information: start, end, color, y position
    # into a plot dic.
    
    plot_dic = {'start': [],
                'end': [],
                'color': [],
                'ypos': [],
                'motif_number': []}
    for event in meme_dic:
        for region in meme_dic[event]:
            start = meme_dic[event][region]['motif_relative_start'][0]
            end = meme_dic[event][region]['motif_relative_end'][0]
            motif_number = meme_dic[event][region]['motif_number'][0]
            # offset start and end depending on region
            start += plot_settings_dic[region]['offset']
            end += plot_settings_dic[region]['offset']
            ypos = event_count/10.0
            try:
                # subtract motif number by 1 to get 0-based numbering
                color = plot_settings_dic[region]['color'][motif_number-1]
            except IndexError:
                print region, motif_number
                print 'Ran out of colors, using yellow as default.'
                color = 'yellow'
            for key, value in \
                zip(['start', 'end', 'color', 'ypos', 'motif_number'], 
                    [start, end, color, ypos, '%s:Motif %s'%(region, motif_number)]):
                plot_dic[key].append(value)
            event_count += 1
            
    if raw_locations is False:
        # begin: get lists of starts, colors, labels for density plot
        density_plot_dic = {}
        for motif_number, color, start in \
            zip(plot_dic['motif_number'], plot_dic['color'], plot_dic['start']):
            if motif_number not in density_plot_dic:
                density_plot_dic[motif_number] = {'densitystarts': []}
            density_plot_dic[motif_number]['densitycolors'] = color
            density_plot_dic[motif_number]['densitystarts'].append(start)
        starts_list = []
        labels_list = []
        colors_list = []
        for motif_number in density_plot_dic:
            labels_list.append(motif_number)
            colors_list.append(density_plot_dic[motif_number]['densitycolors'])
            starts_list.append(density_plot_dic[motif_number]['densitystarts'])
        # add number of sites in labels
        labels_with_nsites = []
        motif_labels = ['Motif %s' %n for n in range(1, len(starts_list) + 1)]
        for labellist, startlist in zip(motif_labels, starts_list):
            n_sites = len(startlist)
            labels_with_nsites.append('%s (%s sites)' %(labellist, n_sites))
        for startlist in starts_list:
            print 'Number of guys: %s' %len(startlist)
            print 'Min/Max: %s/%s' %(min(startlist), max(startlist))
        # end: get lists of starts, colors, labels for density plot
        # begin: init figure
        fig = plt.figure()
        ax = fig.add_subplot(111)
        # end: init figure
        # add rectangles representing exons
        rect_height=0.002
        rect_length=10
        rectstarts = [offsets[0], offsets[2], offsets[4]]
        rcolors = ['cyan', 'yellow', 'cyan']
        for start, color in zip(rectstarts, rcolors):
            patch = add_rectangles(start, height=rect_height, length=rect_length, color=color)
            ax.add_patch(patch)
        # draw intron lines connecting exons
        istarts = [offsets[0], offsets[1], offsets[2], offsets[3]]
        iends = [offsets[1] - rect_length, offsets[2] - rect_length, offsets[3] - rect_length, offsets[4] - rect_length]
        for start, end in zip(istarts, iends):
            plt.hlines(y=-rect_height/2., xmin=start, xmax=end, 
                       color='black', linewidths=1.5)
        # draw vertical lines representing break in intron
        breakstarts = [iends[0], istarts[1], iends[2], istarts[3]]
        for bstart in breakstarts:
            plt.vlines(bstart, ymin=-rect_height, ymax=0, 
                       color='black', linewidths=1)
        plot_functions.plot_density(values_lists=starts_list,
                                    labels_lists=labels_with_nsites,
                                    colors_list=colors_list,
                                    mytitle='Intronic distribution of MEME motifs',
                                    xlabel='Genic region',
                                    ylabel='Density',
                                    xmin=0,
                                    xmax=440,
                                    smoothness=0.1,
                                    legend_pos=2,
                                    showplot=False)
        # dont show xaxis
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.show()
    else:
        plot_functions.plot_hline_segments(starts=plot_dic['start'], 
                                           stops=plot_dic['end'], 
                                           ypos=plot_dic['ypos'], 
                                           colors=plot_dic['color'], 
                                           labels=plot_dic['motif_number'])
    
    
    
if __name__ == '__main__':
    main()