'''
Created on 2013-05-22

@author: jyeung

If there is an AS event: Find the exon/junction with the most variants
If no AS event: use the gene expression
'''

import os
import sys
from utilities import set_directories, find_spliced_events


# Set directories
_cur_dir, \
_proj_dir, \
_input_dir, \
_output_dir = set_directories.set_directories('input', 'takeda', 'output')

'''
_cur_dir = os.path.dirname(os.path.realpath(__file__))
_proj_dir = os.path.dirname(os.path.dirname(os.path.dirname(_cur_dir)))
_input_dir = os.path.join(_proj_dir, 'inputs', 'takeda')
_output_dir = os.path.join(_proj_dir, 'outputs')
_plot_dir = os.path.join(_output_dir, 'plots')
'''

# Set constants
output_path = os.path.join(_output_dir, 'tables')
gene_expr_mean_colname = 'mean'
firmacount_colname = 'w<0.7_count'
ase_group_count_colname = 'count'
gene_id_colname = 'gene_ID'
gene_symbol_colname = 'gene_symbol'
si_sd_colname = 'SI_SD'
probe_id_colname = 'id'
full_sample_list = ['X946_Urethra', 
                    'X972.2.Penila',
                    'AB352',
                    'C42.RNA',
                    'LN.AI.Luc',
                    'X963.Lpa.JP',
                    'X963.L.LN',
                    'X1005',
                    'X890L',
                    'X890.LN',
                    'X945',
                    'X945L.LN',
                    'X961']

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Splicing data must be provided in the command line.')
        sys.exit()
    
    splice_fname = sys.argv[1]
    splice_path = os.path.join(_input_dir, splice_fname)
    
    splice_dat = find_spliced_events.spliced_data(splice_path)
    
    print splice_dat.fpath
    splice_dat.find_spliced_events(output_path, 
                                   gene_expr_mean_colname,
                                   firmacount_colname,
                                   ase_group_count_colname,
                                   gene_id_colname,
                                   gene_symbol_colname,
                                   si_sd_colname,
                                   probe_id_colname,
                                   full_sample_list)
    
    
    
    
    
    
    
    
    