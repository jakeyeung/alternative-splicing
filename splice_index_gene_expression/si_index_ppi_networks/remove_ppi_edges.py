'''
Created on 2013-07-04

@author: jyeung

Remove PPI edges where two genes's expression is below some user defined threshold.
'''


import os
import csv
import sys
from scipy import stats
from utilities import set_directories


# Set directories
mydirs = set_directories.mydirs('input', 'output')


def read_gene_expression(gene_expression_fullpath, samplenames_list, gene_symbol_colname):
    '''
    Read gene expression file, it should contain as column names:
    ensemblID, id, gene_symbol, samplenames, mean, SD.
    
    We will read this file, obtain its gene_symbol and create a dictionary
    with gene_symbol as 'key', list of expression values (ordered by samplenames_list)
    as 'values'. 
    
    Note this file contains multiple probes per gene, but since the file is gene expression
    each probe will have the same value within a gene.
    '''
    # gene_symbol_colname = 'gene_symbol'
    # Init dictionary with gene_symbol as key, list of exprs across samples as value
    gene_expression_dic = {}
    
    with open(gene_expression_fullpath, 'rb') as readfile:
        reader = csv.reader(readfile, delimiter='\t')
        # Read the headers
        colnames = reader.next()
        for row in reader:
            gene_symbol = row[colnames.index(gene_symbol_colname)]
            if gene_symbol not in gene_expression_dic:
                '''
                If gene_symbol not in gene_expression_dic, then we must
                add a list of gene expression across samples into the dic.
                '''
                gene_exprs_list = []
                for s in samplenames_list:
                    gene_exprs_list.append(float(row[colnames.index(s)]))
                gene_expression_dic[gene_symbol] = gene_exprs_list
            else:
                '''
                If gene_symbol already in gene_expression_dic, it means
                we have already, so skip.
                '''
                pass
    return gene_expression_dic


def remove_edges_from_ppi(ppi_fullpath, gene_exprs_dic, 
                          frac_threshold, output_fullpath,
                          frac_samples_to_pass_threshold=1, 
                          use_gene_exprs_threshold=False,
                          normalize_exprs=False,
                          header=False, write_to_file=True):
    '''
    Read each line in ppi fullpath, then check if the two genes is above
    the frac_threshold (to be calculated).
    If both are below frac_threshold,
    do not write the gene pair into output file.
    Else write gene pair into output file. 
    
    Header by default is False, if it contains header, switch to True to skip 
    first row. 
    '''
    
    # BEGIN: Calculate threshold gene expression from frac_threshold.
    # Create long list of gene expression values.
    gene_exprs_list = sorted([i for sublist in gene_exprs_dic.values() \
                              for i in sublist])
    if normalize_exprs==True:
        gene_exprs_list = stats.zscore(gene_exprs_list)
    elif normalize_exprs==False:
        pass
    else:
        sys.exit('normalize_exprs must be True or False, %s detected.' \
                 %normalize_exprs)
    if use_gene_exprs_threshold==False:
        gene_exprs_threshold = gene_exprs_list[int(frac_threshold * len(gene_exprs_list))]
    else:
        gene_exprs_threshold = use_gene_exprs_threshold
    print('log2 gene exprs threshold to pass = %.3f' %gene_exprs_threshold)
    
    # END: Calculate threshold gene exprs from frac_threshold.
    
    # BEGIN: Read PPI network, check if two genes are above threshold, write 
    # to file if above, do not write otherwise.
    reject_count = 0
    row_count = 0
    interactions_not_found = 0
    with open(ppi_fullpath, 'rb') as readfile, \
         open(output_fullpath, 'wb') as writefile:
        reader = csv.reader(readfile, delimiter='\t')
        writer = csv.writer(writefile, delimiter='\t')
        if header == True:
            reader.next()
        elif header == False:
            pass
        else:
            sys.exit('header must be True or False.')
        for row in reader:
            row_count += 1
            '''
            For each row, assumes gene1 is in first column, 
            gene2 is in second column.
            '''
            gene1_name = row[0]
            gene2_name = row[1]
            try:
                gene1_exprs_list = sorted(gene_exprs_dic[gene1_name], reverse=True)
                gene2_exprs_list = sorted(gene_exprs_dic[gene2_name], reverse=True)
                
                '''
                # Define gene1/2 exprs as the maximum value in a gene exprs list that has been
                # filtered to only include a fraction of the samples.
                # This means if frac_samples_to_pass_threshold==0.5, half the samples must be 
                # above the threshold in order to be considered included in the tissue. 
                '''
                gene1_exprs = min(gene1_exprs_list[0:int(frac_samples_to_pass_threshold*len(gene1_exprs_list))])
                gene2_exprs = min(gene2_exprs_list[0:int(frac_samples_to_pass_threshold*len(gene2_exprs_list))])
                if gene1_exprs and gene2_exprs >= gene_exprs_threshold:
                    if write_to_file==True:
                        writer.writerow([gene1_name, gene2_name, 1])
                    elif write_to_file==False:
                        pass
                    else:
                        sys.exit('Writefile must be True or False, %s found'\
                                  %writefile)
                else:
                    if write_to_file==True:
                        writer.writerow([gene1_name, gene2_name, 0])
                    elif write_to_file==False:
                        pass
                    else:
                        sys.exit('Writefile must be True or False, %s found'\
                                  %writefile)
                    reject_count += 1
            except KeyError:
                if gene1_name not in gene_exprs_dic:
                    # print('%s not found, skipping.' %gene1_name)
                    interactions_not_found += 1
                elif gene2_name not in gene_exprs_dic:
                    # print('%s not found, skipping.' %gene2_name)
                    interactions_not_found += 1
                else:
                    sys.exit('%s and %s are strange, human check required.' %(gene1_name, gene2_name))
    # END: Read PPI network...
    print('%s interactions removed out of %s' %(reject_count, row_count-interactions_not_found))
    print('%s genes in PPI not found in cohort out of %s.' \
          %(interactions_not_found, row_count))
    return {'row_count': row_count, 'unmapped_genes': interactions_not_found, 
            'edges_removed': reject_count, 'mapped_genes':row_count-interactions_not_found,
            'gene_exprs_threshold': gene_exprs_threshold}
 

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('Gene expression file, ppi_network, CSV samplenames, '\
              'fractional threshold, and output filename must be indicated in commandline.'\
              'gene expression file and ppi network must be relative to input folder.')
        sys.exit()
    gene_expression_filename = sys.argv[1]
    ppi_filename = sys.argv[2]
    samplenames_csv = sys.argv[3]
    gene_symbol_colname = sys.argv[4]
    try:
        frac_threshold = float(sys.argv[5])
    except ValueError:
        sys.exit('Could not convert %s to float.' %frac_threshold)
    frac_samples_to_pass_threshold = float(sys.argv[6])
    use_gene_exprs_threshold = float(sys.argv[7])
    output_filename = sys.argv[8]
    
    gene_expression_fullpath = os.path.join(mydirs.inputdir, gene_expression_filename)
    ppi_fullpath = os.path.join(mydirs.inputdir, ppi_filename)
    output_fullpath = os.path.join(mydirs.outputdir, output_filename)
    samplenames_list = samplenames_csv.split(',')
    
    # Find threshold by reading gene expression, and saying anything
    # below frac_threshold is too low to be considered expressed.
    gene_exprs_dic = read_gene_expression(gene_expression_fullpath, 
                                          samplenames_list, gene_symbol_colname)
    
    summary_dic = remove_edges_from_ppi(ppi_fullpath, gene_exprs_dic, 
                                        frac_threshold,
                                        output_fullpath,
                                        frac_samples_to_pass_threshold=frac_samples_to_pass_threshold,
                                        use_gene_exprs_threshold=use_gene_exprs_threshold,
                                        normalize_exprs=True, 
                                        header=False, write_to_file=True)
    print summary_dic
    
    
    