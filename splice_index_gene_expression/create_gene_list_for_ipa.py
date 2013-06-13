'''
Created on 2013-06-13

@author: jyeung
create textfile containing 3 columns.
Column 1:
    List of PPI genes
Column 2:
    X or '' depending on if gene showed up in optdis
Column 3:
    log2 fold change between NEPC and PC.
'''


import sys
import csv
import os
from utilities import set_directories


# set constants
inputdir = 'input'
outputdir = 'output'
mydirs = set_directories.mydirs(inputdir, outputdir)
output_colnames = ['background_genes', 'in_list_identifier', 'fold_change']


def get_gene_list_from_multicolumn_file(file_fullpath, gene_columns, header=False):
    '''
    Read PPI file, which contains two columns of genes.
    Obtain list of genes (unique genes only) and return it.
    
    gene_columns specifies which column index to search to obtain the ppi genes. 
    '''
    ppi_genes = []
    with open(file_fullpath, 'rb') as ppi_file:
        ppi_reader = csv.reader(ppi_file, delimiter='\t')
        if header == True:
            ppi_reader.next()
        else:
            pass
        for row in ppi_reader:
            for i in gene_columns:
                ppi_genes.append(row[i])
    # Get unique genes only, return as a list.
    return list(set(ppi_genes))

def match_genes_to_gene_list(genes_to_match, full_gene_list):
    '''
    Returns a list of string of either '' or 'X', where 'X' indicates
    the gene matches 
    '''
    output_list = []    # Either '' or 'X', X means gene matches.
    for bckgrnd_gene in full_gene_list:
        if bckgrnd_gene in genes_to_match:
            output_list.append('X')
        else:
            output_list.append('')
    return output_list

def calculate_fold_change_in_optdis_results(optdis_input_fullpath, group1_genes, group2_genes):
    '''
    Calculates fold change from PPI genes. Requires optdis input file.
    
    Output is a dictionary with genes as keyname, fold change between (group1 / group2)
    as values. 
    '''
    output_dict = {}
    with open(optdis_input_fullpath, 'rb') as optdis_file:
        optdis_reader = csv.reader(optdis_file, delimiter='\t')
        colheaders = optdis_reader.next()    # File contains headers always
        optdis_reader.next()    # Skip one row containing class names. 
        for row in optdis_reader:
            genename = row[0]    # First column
            group1_exprs_list = []
            group2_exprs_list = []
            for g1_sample in group1_genes:
                group1_exprs_list.append(float(row[colheaders.index(g1_sample)]))
            for g2_sample in group2_genes:
                group2_exprs_list.append(float(row[colheaders.index(g2_sample)]))
            group1_mean_exprs = float(sum(group1_exprs_list)) / len(group1_exprs_list)
            group2_mean_exprs = float(sum(group2_exprs_list)) / len(group2_exprs_list)
            # Append to dictionary as fold change between two groups
            output_dict[genename] = float(group1_mean_exprs) - float(group2_mean_exprs)
    return output_dict

def match_fold_change_dict_to_gene_list(gene_fold_change_dict, gene_list):
    fold_change_list = []
    for gene in gene_list:
        try:
            fold_change_list.append(gene_fold_change_dict[gene])
        except KeyError:
            print('%s not found in dict' %gene)
            fold_change_list.append('NA')
    return fold_change_list
    
def write_ppi_genes_optdis_genes_fold_change(ppi_genes, optdis_genes_in_ppi, 
                                             fold_change_list, 
                                             output_colnames,
                                             output_fullpath):
    with open(output_fullpath, 'wb') as output_file:
        output_writer = csv.writer(output_file, delimiter='\t')
        # Write header
        output_writer.writerow(output_colnames)
        for gene, in_list, fc in zip(ppi_genes, optdis_genes_in_ppi, fold_change_list):
            output_writer.writerow([gene, in_list, fc])
    return None

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print('PPI network, optdis output (appended with gene expression)'\
              ' appended with gene/probe expression, NEPC samplenames'\
              ' (CSV format), and PC samplenames (CSV format) must be '\
              ' specified in command'\
              '\nPath names are relative to output folder.')
        sys.exit()
    ppi_partialpath = sys.argv[1]
    optdis_output_partialpath = sys.argv[2]
    optdis_input_partialpath = sys.argv[3]
    output_partialpath = sys.argv[4]
    nepc_genes_csv = sys.argv[5]
    pc_genes_csv = sys.argv[6]
    
    nepc_gene_list = nepc_genes_csv.split(',')
    pc_gene_list = pc_genes_csv.split(',')
    
    ppi_fullpath = os.path.join(mydirs.outputdir, ppi_partialpath)
    optdis_output_fullpath = os.path.join(mydirs.outputdir, optdis_output_partialpath)
    optdis_input_fullpath = os.path.join(mydirs.outputdir, optdis_input_partialpath)
    output_fullpath = os.path.join(mydirs.outputdir, output_partialpath)

    ppi_genes = get_gene_list_from_multicolumn_file(ppi_fullpath, [0, 1], header=False)
    
    optdis_genes = get_gene_list_from_multicolumn_file(optdis_output_fullpath, [0], header=True)
    
    optdis_genes_in_ppi = match_genes_to_gene_list(optdis_genes, ppi_genes)
    
    gene_fold_change_dict = calculate_fold_change_in_optdis_results(optdis_input_fullpath, nepc_gene_list, pc_gene_list)
    
    fold_change_list = match_fold_change_dict_to_gene_list(gene_fold_change_dict, ppi_genes)
    
    write_ppi_genes_optdis_genes_fold_change(ppi_genes, optdis_genes_in_ppi, fold_change_list, output_colnames, output_fullpath)
    
    
    
    
    