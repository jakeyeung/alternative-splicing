'''
Created on 2013-08-30

@author: jyeung
Given a list of gene names (AS events), load a PPI network and calculate
the gene's neighborhood expression.
'''


import csv
import sys
from utilities import igraph_utilities
from miso_scripts import group_miso_utils


def load_genes(filepath, index_list):
    '''
    Given a text file containing only 
    gene names (one gene name per row),
    load the genes to memory.
    '''
    gene_list = []
    with open(filepath, 'rb') as myfile:
        myreader = csv.reader(myfile, delimiter='\t')
        for row in myreader:
            for i in index_list:
                gene_list.append(row[i])
    return gene_list

def load_gene_exprs(gene_exprs_file, samples):
    '''
    Expected as gene symbol colname: 'gene'
    and Sample names (e.g. X97_T, note the X prefix.)
    
    Output dic format:
    genesymbol as key,
        inside is a list of gene expression values, ordered
        the same as samples. 
    '''
    # Define column name constant.
    gene_str = 'gene'
    
    gene_exprs_dic = {}
    with open(gene_exprs_file, 'rb') as myfile:
        myreader = csv.reader(myfile, delimiter='\t')
        colnames = myreader.next()    # Important for getting indices later.
        for row in myreader:
            gene = row[colnames.index(gene_str)]
            # we will not allow duplicate gene entries, so check if in dic.
            if gene not in gene_exprs_dic:
                gene_exprs_dic[gene] = []    # initialize empty list
                for s in samples:
                    gene_exprs_dic[gene].append\
                        ((s, float(row[colnames.index(s)])))
            else:
                print('Warning, duplicate on gene %s, skipping...' %gene)
    return gene_exprs_dic

def get_gene_exprs_two_groups(gene_exprs_dic,
                              g1_samples_list,
                              g2_samples_list):
    '''
    gene_exprs_dic expected to contain
    {gene_symbol: [(s1, gene_exprs_samp1), ... (s30, samp30)]}
    We want instead,
    {gene_symbol: (group1_mean_exprs, group2_mean_exprs)}
    '''
    for gene in gene_exprs_dic:    # iterate over keys
        g1_exprs = []
        g2_exprs = []
        all_exprs = [i[1] for i in gene_exprs_dic[gene]]
        all_samples = [i[0] for i in gene_exprs_dic[gene]]
        for i in range(0, len(all_samples)):
            if all_samples[i] in g1_samples_list:
                g1_exprs.append(float(all_exprs[i]))
            elif all_samples[i] in g2_samples_list:
                g2_exprs.append(float(all_exprs[i]))
        mean_g1 = float(sum(g1_exprs)) / len(g1_exprs)
        mean_g2 = float(sum(g2_exprs)) / len(g2_exprs)
        # Overwrite with simpler values, just a tuple.
        gene_exprs_dic[gene] = (mean_g1, mean_g2)
    return gene_exprs_dic

def calc_abs_exprs_diff(gene, gene_exprs_dic):
    '''
    Given a gene, get their absolute expression difference between
    the two groups.
    
    Note, the specified gene may not be in gene_exprs_dic, skip if so.
    '''
    try:
        abs_diff = abs(gene_exprs_dic[gene][0] - gene_exprs_dic[gene][1]) 
    except KeyError:
        abs_diff = None
    return abs_diff

def main():
    if len(sys.argv) < 3:
        print('Candidate gene file, gene exprs file and PPI network must be' \
              ' specified in command line.')
        sys.exit()
    gene_names_file = sys.argv[1]
    gene_exprs_file = sys.argv[2]
    ppi_file = sys.argv[3]
    pc_samples_file = sys.argv[4]
    nepc_samples_file = sys.argv[5]
    output_file = sys.argv[6]
    
    # Define write colname constants
    as_gene_str = 'as_gene'
    neighbors_str = 'neighbors'
    n_neighbors_str = 'number_of_neighbors'
    sum_abs_diff_str = 'absolute_gene_exprs_diff_score'
    sum_abs_diff_norm_str = 'absolute_gene_exprs_diff_score(degree-normalized)'
    as_fold_change_str = 'as_gene_fold_change'
    colnames = [as_gene_str, neighbors_str, n_neighbors_str, sum_abs_diff_str, 
                sum_abs_diff_norm_str, as_fold_change_str]
    
    # Load gene names, genes of interest, such as AS genes.
    # Index list 0 because it is first row only.
    gene_list = load_genes(gene_names_file, index_list=[0])
    gene_list_length_total = len(gene_list)
    
    # Index list 0, 1 because we want both rows.
    ppi_genes = load_genes(ppi_file, index_list=[0, 1])
    
    # Get intersection of gene_list and ppi_genes, we unforunately cannot
    # do calculations on genes that do not overlap with ppi genes.
    gene_list = list(set(gene_list).intersection(set(ppi_genes)))
    gene_list_length_subset = len(gene_list)
    print('%s of %s genes removed because they do not map to PPI genes.' \
            %(gene_list_length_subset, gene_list_length_total))
    
    # Load PPI network as graph object graph
    graph = igraph_utilities.igraph_object(ppi_file)
    
    # Load samples to two groups.
    pc_samples_list = \
        group_miso_utils.get_sample_names_from_file(pc_samples_file)
    nepc_samples_list = \
        group_miso_utils.get_sample_names_from_file(nepc_samples_file)
        
    # Load gene expression into dictionary
    all_samples = pc_samples_list + nepc_samples_list
    gene_exprs_dic = load_gene_exprs(gene_exprs_file, all_samples)
    
    '''
    # Currently ,the dic contains list of gene exprs for all samples
    # let's shorten it by getting mean exprs of the two groups.
    # Modify gene_exprs_dic so each value contains a list of only two 
    # values, mean exprs of pc and mean exprs of nepc (respectively).
    '''
    gene_exprs_dic = get_gene_exprs_two_groups(gene_exprs_dic, 
                                               pc_samples_list, 
                                               nepc_samples_list)
    # Find and calculate neighbor gene expression around specified genes.
    # Write to outputs to file.
    with open(output_file, 'wb') as myfile:
        writer = csv.writer(myfile, delimiter='\t')
        # Write header
        writer.writerow(colnames)
        
        for gene in gene_list:
            # Get fold change of gene of interest
            as_fold_change = calc_abs_exprs_diff(gene, gene_exprs_dic)
            
            # Get its neighbors and its node degree.
            neighbor_genes = graph.get_neighbors(gene, include_self=False)
            # Calculate expression difference
            neighbor_genes_with_exprs = []
            abs_exprs_diff_list = []
            neighbor_count = 0    # init counter
            for neighbor_g in neighbor_genes:
                abs_exprs_diff = calc_abs_exprs_diff(neighbor_g, 
                                                gene_exprs_dic)
                if abs_exprs_diff != None:
                    '''
                    If None, it means gene exprs could not be found for
                    that gene. 
                    Otherwise, we will append to list, add neighbor_count
                    and then expand list of neighbors that contain
                    gene exprs (useful for normalizing and writing to file)
                    '''
                    abs_exprs_diff_list.append(abs_exprs_diff)
                    neighbor_count += 1
                    neighbor_genes_with_exprs.append(neighbor_g)              
            # Sum the values up.
            abs_exprs_sum = sum(abs_exprs_diff_list)
            try:
                abs_exprs_sum_norm = abs_exprs_sum / len(abs_exprs_diff_list)
            except ZeroDivisionError:
                # Empty list, then call it NA
                abs_exprs_sum_norm = 'NA'
            # Write to file, convert any lists to a csv val
            myrow = [gene, ','.join(neighbor_genes_with_exprs),
                     neighbor_count,
                     abs_exprs_sum, abs_exprs_sum_norm, as_fold_change]
            writer.writerow(myrow)
    print('Output file saved to %s' %output_file)
            
if __name__ == '__main__':
    main()