'''
Created on 2013-05-22

@author: jyeung
'''


import csv
import sys
import dict_tools


class spliced_data:
    '''
    Class representing spliced data
    '''
    
    def __init__(self, full_path_name):
        '''
        Initialize
        '''
        self.fpath = full_path_name
    
    def find_spliced_events(self, output_path, 
                            gene_expr_mean_colname, 
                            firmacount_colname, 
                            ase_group_count_colname,
                            gene_id_colname,
                            gene_symbol_colname,
                            si_sd_colname,
                            probe_id_colname,
                            full_sample_list):
        '''
        gene_expr_mean_colname = 'mean'
        w_count_colname = 'w<0.7_count'
        ASE_group_count = 'count'
        
        Beware of duplicated column names!
        
        Read splicing data, check for highly confident ASE. If ASE, take 
        exon/junction of highest variation, else take mean gene expression. 
        '''
        with open(self.fpath, 'rb') as splice_file:
            # Open file
            splice_reader = csv.reader(splice_file, delimiter='\t')
            # Define columns
            read_colnames = splice_reader.next()
            mean_gexp_index = read_colnames.index(gene_expr_mean_colname)
            firmacount_index = read_colnames.index(firmacount_colname)
            ase_group_count_index = read_colnames.index(ase_group_count_colname)
            gene_id_index = read_colnames.index(gene_id_colname)
            gene_symbol_index = read_colnames.index(gene_symbol_colname)
            si_sd_index = read_colnames.index(si_sd_colname)
            probe_id_index = read_colnames.index(probe_id_colname)
            
            # Define sample columns as the indexes from 
            # si_sd_index - len(sample_list) to sd_sd_index - 1
            sample_si_indices_start = si_sd_index - len(full_sample_list)
            sample_si_index_list = range(sample_si_indices_start, si_sd_index)
            
            # Define gene expression columns as the indexes from
            # mean - len(sample_list) to mean - 1
            gene_exprs_indices_start = mean_gexp_index - len(full_sample_list)
            gene_exprs_index_list = range(gene_exprs_indices_start, 
                                          mean_gexp_index)
            # Check the indexes are correct by matching full_sample_list with 
            # column name.
            sidata_sample_colnames = [read_colnames[i] for i in sample_si_index_list]
            genedata_sample_colnames = [read_colnames[i] for i in gene_exprs_index_list]
            
            # Check the sidata_sample_colnames == full_sample_list
            if sidata_sample_colnames == full_sample_list:
                if genedata_sample_colnames == full_sample_list:
                    print('Splicing data and gExp sample colnames found...')
                else:
                    print('gExp sample colnames do not match user input sample list')
                    sys.exit('Exiting...')
            else:
                print('Data sample colnames do not match user input sample list')
                sys.exit('Exiting...')
            # BEGIN: Create key names with an empty list for dictionaries
            # Initialize dicts to store final data (one gene, one number)
            gene_info = {}
            # Initialize dicts to store exon-level data (within one gene)
            # Values from exon_info will be used to enlargen the list dictionaries
            # of genes_with_as and genes_without_as.
            exon_info = {}
            # Add suffix _si or _mean_gene_exprs to sidata_sample_colnames
            sample_si_names = ['%s_si' %i for i in full_sample_list]
            sample_mean_exprs_names = ['%s_mean_exprs' \
                                            %i for i in full_sample_list]
            # Create list of key names to be added to dictionary with empty list
            gene_keynames = ['gene_id', 'gene_symbol', 
                             'probe'] + full_sample_list
            # non_as_keynames = ['gene_id', 'gene_symbol'] + sample_mean_exprs_names
            exon_info_keynames = ['gene_id', 
                                  'gene_symbol',
                                  'probe',
                                  'si_sd'] + sample_si_names + sample_mean_exprs_names
            # Create key:emptylist keyval pairs
            for gkey in gene_keynames:
                gene_info = dict_tools.init_key_emptylist(gene_info, 
                                                          gkey)
            for exon_key in exon_info_keynames:
                exon_info = dict_tools.init_key_emptylist(exon_info, 
                                                          exon_key)
            # END: Create key names with an empty list for gene-wide dictionary
            
            # BEGIN: Create dict for checking if the exon/junctions in the gene
            # exhibit any AS events. Values here will not be transferred 
            # to genes_with_as or genes_without_as dicts.. 
            exon_filters = {}
            # Create key:emptylist keyval pairs
            filter_keynames = ['ASE_group_count', 'firma_count', 'mean_gene_exprs']
            for filter_key in filter_keynames:
                exon_filters = dict_tools.init_key_emptylist(exon_filters, filter_key)
                
            # END: Create dict for checking if the exon/junctions in the gene
            # exhibit any AS events.
            
            # BEGIN: Iterate each row to separate AS events with non AS events
            # Initialize some values that will be used
            rowcount = 0
            prev_gene = 0
                
            for row in splice_reader:
                '''
                Iterate rows to collect exon information.
                Before each row collection, check that prev_gene == curr_gene
                if they are not equal, then we've collected all the 
                exon information for that gene. We go into an if loop 
                to check if this gene has ASE. 
                If ASE occurred, then we take most varied probe expression.
                If no ASE, we take gene expression. 
                '''
                curr_gene = row[gene_id_index]

                if curr_gene != prev_gene and rowcount > 0:
                    print('%s rows collected, analyzing info...' %rowcount)
                    '''
                    If curr_gene != prev_gene, we've collected all exon info
                    for this gene. Time to analyze results from our 
                    exon-level dicts and put it intoour gene-level dicts.
                    
                    First we check if splicing occurred by seeing if:
                    mean_gexprs >= 3 for all exons
                    ASEgroupcount >= 2 for any exon/junction
                    firmaweightcount >= 2 for any exon/junction.
                    
                    If splicing occurred, find most varied probe across
                    all samples (max si_sd).
                    Using most varied probe, fill out geneID, geneSymbol,
                    probe, and probe expression of each sample. 
                    
                    If no splicing occurred, then fill out geneID, geneSymbol,
                    set probe to all_probes, and expression of entire gene of 
                    each sample. 
                    '''
                    if min(exon_filters['mean_gene_exprs']) >= 3\
                    and min(exon_filters['ASE_group_count']) >= 2\
                    and min(exon_filters['firma_count']) >= 2:
                        # Find index at max variance
                        maxprobe_index = exon_info['si_sd'].index(max(exon_info['si_sd']))
                        print maxprobe_index
                        gene_info['gene_id'].append(exon_info['gene_id'][maxprobe_index])
                        gene_info['gene_symbol'].append(exon_info['gene_symbol'][maxprobe_index])
                        gene_info['probe'].append(exon_info['probe'][maxprobe_index])
                        for i in range(len(full_sample_list)):
                            gene_info[full_sample_list[i]].append(exon_info[sample_si_names[i]][maxprobe_index])
                    else:
                        # Find index at max variance
                        maxprobe_index = max(exon_info['si_sd'])
                        gene_info['gene_id'].append(exon_info['gene_id'][maxprobe_index])
                        gene_info['gene_symbol'].append(exon_info['gene_symbol'][maxprobe_index])
                        gene_info['probe'].append('all_probes')
                        for i in range(full_sample_list):
                            gene_info[full_sample_list[i]].append(exon_info[sample_mean_exprs_names[i]])
                            
                    # Reinitialize exon_info and exon_filters
                    for key in exon_filters:
                        exon_filters[key] = []
                    for key in exon_info:
                        exon_info[key] = []
                    print('One gene done, rowcount = %s' %rowcount)
                    break
                # Extend lists in dictionary
                exon_info['gene_id'].append(row[gene_id_index])
                exon_info['gene_symbol'].append(row[gene_symbol_index])
                exon_info['probe'].append(row[probe_id_index])
                exon_info['si_sd'].append(row[si_sd_index])
                # BEGIN: Add append SI into exon_info
                samp_count = 0
                for i in sample_si_index_list:
                    exon_info[sample_si_names[samp_count]].append(row[i])
                    samp_count += 1
                # END: Add append SI into exon_info
                # BEGIN: Add append gExp data into exon_info
                samp_count = 0
                for i in gene_exprs_index_list:
                    exon_info[sample_mean_exprs_names[samp_count]].append(row[i])
                    samp_count += 1
                # END: Add append gExp into exon_info
                # BEGIN: Append filter info
                exon_filters['ASE_group_count'].append(row[ase_group_count_index])
                exon_filters['firma_count'].append(row[firmacount_index])
                exon_filters['mean_gene_exprs'].append(row[mean_gexp_index])
                prev_gene = curr_gene
                rowcount += 1
            print gene_info
            # END: Iterate each row to separate AS events with non AS events
                
            
        
        