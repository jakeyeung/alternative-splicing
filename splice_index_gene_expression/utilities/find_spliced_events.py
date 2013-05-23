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
    
    def find_spliced_events_si(self, output_path, 
                               gene_expr_mean_colname, 
                               firmacount_colname, 
                               ase_group_count_colname,
                               gene_id_colname,
                               gene_symbol_colname,
                               si_sd_colname,
                               probe_id_colname,
                               full_sample_list,
                               printlong=False):
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
                             'probe',
                             'probe_or_gene'] + full_sample_list
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
            gene_count = 0
            ase_count = 0
            non_ase_count = 0
                
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
                    '''
                    print exon_filters['mean_gene_exprs']
                    print exon_filters['ASE_group_count']
                    print exon_filters['firma_count']
                    raw_input('Press enter to continue: ')
                    '''
                    if min(exon_filters['mean_gene_exprs']) >= 3\
                    and max(exon_filters['ASE_group_count']) >= 2\
                    and max(exon_filters['firma_count']) >= 2:
                        # Find index at max variance
                        maxprobe_index = exon_info['si_sd'].index(max(exon_info['si_sd']))
                        # print maxprobe_index
                        gene_info['gene_id'].append(exon_info['gene_id'][maxprobe_index])
                        gene_info['gene_symbol'].append(exon_info['gene_symbol'][maxprobe_index])
                        gene_info['probe'].append(exon_info['probe'][maxprobe_index])
                        gene_info['probe_or_gene'].append('probe')
                        # Add SI information to each sample. 
                        for i in range(len(full_sample_list)):
                            gene_info[full_sample_list[i]].append(exon_info[sample_si_names[i]][maxprobe_index])
                        ase_count += 1
                    else:
                        # Find index at max variance
                        maxprobe_index = exon_info['si_sd'].index(max(exon_info['si_sd']))
                        gene_info['gene_id'].append(exon_info['gene_id'][maxprobe_index])
                        gene_info['gene_symbol'].append(exon_info['gene_symbol'][maxprobe_index])
                        gene_info['probe'].append('all_probes')
                        gene_info['probe_or_gene'].append('gene')
                        # Add mean gene expression info to each sample.
                        for i in range(len(full_sample_list)):
                            gene_info[full_sample_list[i]].append(exon_info[sample_mean_exprs_names[i]][maxprobe_index])
                        non_ase_count += 1
                    # Reinitialize exon_info and exon_filters
                    for key in exon_filters:
                        exon_filters[key] = []
                    for key in exon_info:
                        exon_info[key] = []
                    # print('%s gene done, rowcount = %s' %(curr_gene, rowcount))
                    gene_count += 1
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
                # print row[ase_group_count_index]
                # print int(row[ase_group_count_index])
                exon_filters['ASE_group_count'].append(int(row[ase_group_count_index]))
                exon_filters['firma_count'].append(int(row[firmacount_index]))
                exon_filters['mean_gene_exprs'].append(float(row[mean_gexp_index]))
                prev_gene = curr_gene
                rowcount += 1
            print('All genes done. Total rowcount: %s ' %rowcount)
            # END: Iterate each row to separate AS events with non AS events
            # BEGIN: Write to file.
        
        # Decide whether to print output as LONG or SHORT.
        # Short means one gene per row, with many sample columns.
        # Long means the sampleIDs are now in the rows, not columns.
        if printlong == False:
            # Code for writing a table where column names contain sampleIDs.
            with open(output_path, 'wb') as writefile:
                print gene_keynames
                print gene_info.keys()
                splice_writer = csv.DictWriter(writefile, delimiter='\t', 
                                               fieldnames=gene_keynames)
                splice_writer.writeheader()
                for i in xrange(gene_count):
                    rowdict = {}
                    for key, val in gene_info.iteritems():
                        rowdict[key] = val[i]
                    splice_writer.writerow(rowdict)
        
        elif printlong == True:
            # Code for writing long-form. SampleID is a header.   
            with open(output_path, 'wb') as writefile:
                print gene_keynames
                print gene_info.keys()
                fieldnames = ['gene_id', 'gene_symbol', 'probe', 'probe_or_gene', 
                              'sampleID', 'probe_or_geneExp']
                splice_writer = csv.DictWriter(writefile, delimiter='\t', 
                                               fieldnames=fieldnames)
                splice_writer.writeheader()
                for i in xrange(gene_count):
                    rowdict = {}
                    gene_mainkeys = ['gene_id', 'gene_symbol', 'probe',
                                    'probe_or_gene']    # Constant per sample per gene
                    for mainkey in gene_mainkeys:
                        rowdict[mainkey] = gene_info[mainkey][i]
                    for sampkey in full_sample_list:
                        rowdict['sampleID'] = sampkey
                        rowdict['probe_or_geneExp'] = gene_info[sampkey][i]
                        splice_writer.writerow(rowdict)
        else:
            sys.exit('printlong must be either True or False, exitnig...')
        return gene_info
    
    def find_spliced_events_probes(self, output_path, 
                                   gene_expr_mean_colname, 
                                   firmacount_colname, 
                                   ase_group_count_colname,
                                   gene_id_colname,
                                   gene_symbol_colname,
                                   si_sd_colname,
                                   probe_sd_colname,
                                   probe_id_colname,
                                   full_sample_list,
                                   printlong=False):
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
            probe_sd_index = read_colnames.index(probe_sd_colname)
            probe_id_index = read_colnames.index(probe_id_colname)
            
            # Define si sample columns as the indexes from 
            # si_sd_index - len(sample_list) to sd_sd_index - 1
            sample_si_indices_start = si_sd_index - len(full_sample_list)
            sample_si_index_list = range(sample_si_indices_start, si_sd_index)
            
            # Define gene expression sample columns as the indexes from
            # mean - len(sample_list) to mean - 1
            gene_exprs_indices_start = mean_gexp_index - len(full_sample_list)
            gene_exprs_index_list = range(gene_exprs_indices_start, 
                                          mean_gexp_index)
            # Define probe expression sample columns as indexes from
            # exp_sd_index - len(full_sample_list. 
            probe_exprs_indices_start = probe_sd_index - len(full_sample_list)
            probe_exprs_index_list = range(probe_exprs_indices_start, 
                                           probe_sd_index)
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
            sample_probe_exprs_names = ['%s_probe' %i for i in full_sample_list]
            # Create list of key names to be added to dictionary with empty list
            gene_keynames = ['gene_id', 'gene_symbol', 
                             'probe',
                             'probe_or_gene'] + full_sample_list
            # non_as_keynames = ['gene_id', 'gene_symbol'] + sample_mean_exprs_names
            exon_info_keynames = ['gene_id', 
                                  'gene_symbol',
                                  'probe',
                                  'si_sd'] + sample_si_names + \
                                  sample_mean_exprs_names + \
                                  sample_probe_exprs_names
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
            gene_count = 0
            ase_count = 0
            non_ase_count = 0
                
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
                    '''
                    print exon_filters['mean_gene_exprs']
                    print exon_filters['ASE_group_count']
                    print exon_filters['firma_count']
                    raw_input('Press enter to continue: ')
                    '''
                    if min(exon_filters['mean_gene_exprs']) >= 3\
                    and max(exon_filters['ASE_group_count']) >= 2\
                    and max(exon_filters['firma_count']) >= 2:
                        # Find index at max variance
                        maxprobe_index = exon_info['si_sd'].index(max(exon_info['si_sd']))
                        # print maxprobe_index
                        gene_info['gene_id'].append(exon_info['gene_id'][maxprobe_index])
                        gene_info['gene_symbol'].append(exon_info['gene_symbol'][maxprobe_index])
                        gene_info['probe'].append(exon_info['probe'][maxprobe_index])
                        gene_info['probe_or_gene'].append('probe')
                        # Add probe expression information to each sample. 
                        for i in range(len(full_sample_list)):
                            gene_info[full_sample_list[i]].append(exon_info[sample_probe_exprs_names[i]][maxprobe_index])
                        ase_count += 1
                    else:
                        # Find index at max variance
                        maxprobe_index = exon_info['si_sd'].index(max(exon_info['si_sd']))
                        gene_info['gene_id'].append(exon_info['gene_id'][maxprobe_index])
                        gene_info['gene_symbol'].append(exon_info['gene_symbol'][maxprobe_index])
                        gene_info['probe'].append('all_probes')
                        gene_info['probe_or_gene'].append('gene')
                        # Add gene expression information to each sample
                        for i in range(len(full_sample_list)):
                            gene_info[full_sample_list[i]].append(exon_info[sample_mean_exprs_names[i]][maxprobe_index])
                        non_ase_count += 1
                    # Reinitialize exon_info and exon_filters
                    for key in exon_filters:
                        exon_filters[key] = []
                    for key in exon_info:
                        exon_info[key] = []
                    # print('%s gene done, rowcount = %s' %(curr_gene, rowcount))
                    gene_count += 1
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
                # BEGIN: Append probe expression into exon_info
                samp_count = 0
                for i in probe_exprs_index_list:
                    exon_info[sample_probe_exprs_names[samp_count]].append(row[i])
                    samp_count += 1
                # END: Append probe expression into exon_info
                # BEGIN: Append filter info
                # print row[ase_group_count_index]
                # print int(row[ase_group_count_index])
                exon_filters['ASE_group_count'].append(int(row[ase_group_count_index]))
                exon_filters['firma_count'].append(int(row[firmacount_index]))
                exon_filters['mean_gene_exprs'].append(float(row[mean_gexp_index]))
                prev_gene = curr_gene
                rowcount += 1
            print('All genes done. Total rowcount: %s ' %rowcount)
            # END: Iterate each row to separate AS events with non AS events
            # BEGIN: Write to file.
        
        # Decide whether to print output as LONG or SHORT.
        # Short means one gene per row, with many sample columns.
        # Long means the sampleIDs are now in the rows, not columns.
        if printlong == False:
            # Code for writing a table where column names contain sampleIDs.
            with open(output_path, 'wb') as writefile:
                print gene_keynames
                print gene_info.keys()
                splice_writer = csv.DictWriter(writefile, delimiter='\t', 
                                               fieldnames=gene_keynames)
                splice_writer.writeheader()
                for i in xrange(gene_count):
                    rowdict = {}
                    for key, val in gene_info.iteritems():
                        rowdict[key] = val[i]
                    splice_writer.writerow(rowdict)
        
        elif printlong == True:
            # Code for writing long-form. SampleID is a header.   
            with open(output_path, 'wb') as writefile:
                print gene_keynames
                print gene_info.keys()
                fieldnames = ['gene_id', 'gene_symbol', 'probe', 'probe_or_gene', 
                              'sampleID', 'probe_or_geneExp']
                splice_writer = csv.DictWriter(writefile, delimiter='\t', 
                                               fieldnames=fieldnames)
                splice_writer.writeheader()
                for i in xrange(gene_count):
                    rowdict = {}
                    gene_mainkeys = ['gene_id', 'gene_symbol', 'probe',
                                    'probe_or_gene']    # Constant per sample per gene
                    for mainkey in gene_mainkeys:
                        rowdict[mainkey] = gene_info[mainkey][i]
                    for sampkey in full_sample_list:
                        rowdict['sampleID'] = sampkey
                        rowdict['probe_or_geneExp'] = gene_info[sampkey][i]
                        splice_writer.writerow(rowdict)
        else:
            sys.exit('printlong must be either True or False, exitnig...')
        return gene_info
            
        
        