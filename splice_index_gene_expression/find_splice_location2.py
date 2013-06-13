'''
Created on 2013-06-12

@author: jyeung
Improvements from find_splice_location.py:
    1) Want to use SI values only to detect whether an AS event happened or not. 
    2) Want to use Probe Exps to select WHICH probes to choose from, the reasoning behind
        this is that SI values is calculated based on the entire sample set, but if we are
        looking at only a subset of these samples, we want to choose probes that are more
        indicative of 'separation' between the two samples.
    3) Read from three separate files to prevent the issue of duplicate column names.
'''


import sys
import os
import numpy as np
from utilities import set_directories, read_write_multiple_files

# Set folder constants
inputdir = 'input'
outputdir = 'output'
mydirs = set_directories.mydirs(inputdir, outputdir)

# Set column name constants
gene_exprs_str = 'gene_exprs'
probe_exprs_str = 'probe_exprs'
si_values_str = 'si_values'
gene_id_colname = 'ensemblID'
probe_id_colname = 'id'
gene_symbol_colname = 'gene_symbol'
probe_mean_colname = 'EXP_mean'
gene_sd_colname = 'SD'
probe_sd_colname = 'EXP_SD'

# Filter criteria column names
gene_mean_colname = 'mean'
firmacount_colname = 'w<0_7_count'
asegroupcount_colname = 'count'

# Samples
sample_list = ['X946_Urethra', 'X972_2_Penila', 'AB352', 'X1005', 'X890L', 'X945', 'X961']


def find_and_choose_spliced_probes(gene_exprs_fullpath, probe_exprs_fullpath, 
                                   si_values_fullpath, output_fullpath, 
                                   gene_exprs_str, probe_exprs_str, si_values_str, 
                                   gene_id_colname,
                                   probe_id_colname,
                                   gene_symbol_colname, 
                                   asegroupcount_colname, 
                                   firmacount_colname, 
                                   gene_mean_colname, 
                                   sample_list):
    '''
    Read gene, probe and si files, find out which genes are alternatively spliced.
    If gene is alternatively spliced, then replace gene expression with probe
    expression that has the highest variance. 
    
    Each row that is read contains row from aech of the three files, so it is contained in a list
    of length 3.
    Therefoer, assume in readrows that first index is gene_exprs file, second index is probe_exprs file
    and third index is si_values file. 
    '''
    # Set constants
    ase_group_count_str = 'ASE_group_count'
    firma_count_str = 'firma_count'
    mean_gene_exprs_str = 'mean_gene_exprs'
    gene_id_str = 'gene_id'
    gene_symbol_str = 'gene_symbol'
    probe_str = 'probe'
    
    criterialist = [ase_group_count_str, firma_count_str, mean_gene_exprs_str]
    read_nicknames = [gene_exprs_str, probe_exprs_str, si_values_str]
    
    sample_list_gene_exprs_appended = ['%s_%s' \
                                       %(s, mean_gene_exprs_str) for s in sample_list]
    sample_list_probe_exprs_appended = ['%s_%s' \
                                        %(s, probe_exprs_str) for s in sample_list]
    
    # Initialize some dictionaries to be used when we iterate through rows.
    filter_criteria = {}
    for criterion in criterialist:
        filter_criteria[criterion] = []

    probe_exprs_dict = {}
    for keyname in [gene_id_str, gene_symbol_str, probe_str] + \
                    sample_list_probe_exprs_appended + \
                    sample_list_gene_exprs_appended:
        probe_exprs_dict[keyname] = []
    
    # Initialize object for reading and writing.
    read_write_obj = read_write_multiple_files.read_write([gene_exprs_fullpath, 
                                                           probe_exprs_fullpath, 
                                                           si_values_fullpath], 
                                                          output_fullpath)
    with read_write_obj:
        # Read header then write header.
        headers_dict = {}
        for name, header in zip(read_nicknames, 
                                read_write_obj.readnext()):
            headers_dict[name] = header
        gene_headers = headers_dict[gene_exprs_str]
        probe_headers = headers_dict[probe_exprs_str]
        si_headers = headers_dict[si_values_str]
        
        # Write header, order is important!
        read_write_obj.writenext([gene_id_str, 
                                 gene_symbol_str, 
                                 probe_str]
                                 + sample_list)
        
        # Iterate through rows, storing information
        rowcount = 0
        prev_gene = 0
        ase_count = 0
        non_ase_count = 0
        while True:
            try:
                rowlist = read_write_obj.readnext()
                rowdict = {}
                for name, rowvals in zip(read_nicknames,
                                        rowlist):
                    rowdict[name] = rowvals
                cur_gene = rowdict[gene_exprs_str][headers_dict[gene_exprs_str].index(gene_symbol_colname)]
            except StopIteration:
                print('%s ASE events, %s non ASE events, breaking...' \
                      %(ase_count, non_ase_count))
                break
            
            if cur_gene != prev_gene and rowcount > 0:
                '''
                If cur_gene does not equal prev_gene, we have read all the probes
                in one gene, so we will gather our information to determine:
                    1) probe of highest variability in that gene
                Then we will write either gene expression (if no ase) or 
                probe expression (if ase detected, use probe of highest variability).
                '''
                if min(filter_criteria[mean_gene_exprs_str]) >= 3\
                and max(filter_criteria[ase_group_count_str]) >= 2\
                and max(filter_criteria[firma_count_str]) >= 2:
                    # Calculate probe expression variance across samples at each probe.
                    probe_var_tuplelist = []
                    
                    for i in range(len(probe_exprs_dict[probe_str])):
                        probe = probe_exprs_dict[probe_str][i]
                        probe_exprs = []
                        for s in sample_list_probe_exprs_appended:
                            probe_exprs.append(probe_exprs_dict[s][i])
                        probe_var = np.std(probe_exprs)
                        probe_var_tuplelist.append((probe_var, i, probe))
                    
                    # Sort list of tuples by probe variance largest to smallest then save
                    # first element in list as variables.
                    probe_var_tuplelist = sorted(probe_var_tuplelist, reverse=True)
                    top_probe_index = probe_var_tuplelist[0][1]
                    
                    # Write highest probe expression to file, and indicate which probe 
                    # was used.
                    writerow = []
                    for keyname in ([gene_id_str, gene_symbol_str, probe_str] + \
                                    sample_list_probe_exprs_appended):
                        writerow.append(probe_exprs_dict[keyname][top_probe_index])
                    read_write_obj.writenext(writerow)
                    ase_count += 1
                else:
                    '''
                    No ASE detected, write only gene expression information.
                    '''
                    writerow = []
                    for keyname in ([gene_id_str, gene_symbol_str, probe_str] + \
                                    sample_list_gene_exprs_appended):
                        if keyname == probe_str:
                            writerow.append('all_probes')
                        else:
                            # Pick any index, doesnt matter they are all same
                            writerow.append(probe_exprs_dict[keyname][0])
                    read_write_obj.writenext(writerow)
                    non_ase_count += 1
                    
                # Reinitialize dicts
                for key in filter_criteria.keys():
                    filter_criteria[key] = []
                for key in probe_exprs_dict.keys():
                    probe_exprs_dict[key] = []
                
            # Extend lists in filter dict
            
            '''
            for filter_dicname, file_colname in zip(criterialist, 
                                                    [asegroupcount_colname, 
                                                     firmacount_colname, 
                                                     gene_mean_colname]):
                filter_criteria[filter_dicname] = \
                rowdict[si_values_str][headers_dict[si_values_str].index(file_colname)]
            '''
            
            '''
            filter_criteria[ase_group_count_str] = \
                rowdict[si_values_str][headers_dict[si_values_str].index(asegroupcount_colname)]
              
            filter_criteria[firma_count_str] = \
                rowdict[si_values_str][headers_dict[si_values_str].index(firmacount_colname)]
                
            filter_criteria[mean_gene_exprs_str] = \
                rowdict[gene_exprs_str][headers_dict[gene_exprs_str].index(gene_mean_colname)]
            '''
            
            # Extend lists filter information dictionary.
            filter_criteria[ase_group_count_str].append(\
                int(rowdict[si_values_str][si_headers.index(asegroupcount_colname)]))
            filter_criteria[firma_count_str].append(\
                int(rowdict[si_values_str][si_headers.index(firmacount_colname)]))
            filter_criteria[mean_gene_exprs_str].append(\
                float(rowdict[gene_exprs_str][gene_headers.index(gene_mean_colname)]))
            
            # Extend lists in probe expression.
            probe_exprs_dict[gene_id_str].append(\
                rowdict[probe_exprs_str][probe_headers.index(gene_id_colname)])
            probe_exprs_dict[gene_symbol_str].append(\
                rowdict[probe_exprs_str][probe_headers.index(gene_symbol_colname)])
            probe_exprs_dict[probe_str].append(\
                rowdict[probe_exprs_str][probe_headers.index(probe_id_colname)])
            # Extend probe and gene expressions for specified samples.
            for s, probe_exprs_key, gene_exprs_key in zip(sample_list, 
                                                          sample_list_probe_exprs_appended, 
                                                          sample_list_gene_exprs_appended):
                probe_exprs_dict[probe_exprs_key].append(\
                    float(rowdict[probe_exprs_str][probe_headers.index(s)]))
                probe_exprs_dict[gene_exprs_key].append(\
                    float(rowdict[gene_exprs_str][gene_headers.index(s)]))
            
            prev_gene = cur_gene
            rowcount += 1
        print('%s rows read, %s rows written, function end.' \
              %(read_write_obj.readrowcount, read_write_obj.writerowcount))
    return None


if __name__ == '__main__':
    if len(sys.argv) < 6:
        print('AS_candidates_gene_exps, AS_candidates_probe_exps, '\
              'AS_candidates_si_values partial paths, Output partial path, '\
              ',comma-separated sample names, and output partial path'\
              ' must be specified in command line.' \
              '\nDatafile paths are relative to the input directory, '\
              'output path is relative to output directory.')
        sys.exit()
    gene_exprs_partialpath = sys.argv[1]
    probe_exprs_partialpath = sys.argv[2]
    si_values_partialpath = sys.argv[3]
    sample_names_string = sys.argv[4]
    output_partialpath = sys.argv[5]
    
    # Fix up user inputs into function inputs.
    gene_exprs_fullpath = os.path.join(mydirs.inputdir, gene_exprs_partialpath)
    probe_exprs_fullpath = os.path.join(mydirs.inputdir, probe_exprs_partialpath)
    si_values_fullpath = os.path.join(mydirs.inputdir, si_values_partialpath)
    sample_names_list = sample_names_string.split(',')
    output_fullpath = os.path.join(mydirs.outputdir, output_partialpath)
    
    # Find when splicing event happens...
    read_write_obj = read_write_multiple_files.read_write([gene_exprs_fullpath, 
                                                           probe_exprs_fullpath, 
                                                           si_values_fullpath], 
                                                          output_fullpath)
    
    find_and_choose_spliced_probes(gene_exprs_fullpath, probe_exprs_fullpath, 
                                   si_values_fullpath, output_fullpath, 
                                   gene_exprs_str, probe_exprs_str, si_values_str, 
                                   gene_id_colname, probe_id_colname, gene_symbol_colname, 
                                   asegroupcount_colname, firmacount_colname, 
                                   gene_mean_colname, sample_list)
    
    
    
    
    
    
    
    
    