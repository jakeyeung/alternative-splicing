'''
Created on 2013-09-27

@author: jyeung
For Raunak's paper:
Find outlier gene expression of NEPC by using
PCa as "null" distribution.
Basically find SD and mean of gene1 in PCa samples,
subtract NEPC sample with PCa mean, ask if it is 
past a certain SD, if so, then it is outlier gene.
'''

import sys
import numpy as np
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
from find_neighborhood_gene_exprs import load_gene_exprs
from utilities import file_utilities

def plot_exprs(exprs, gene, abs_delta_mean, samp_sd):
    '''
    Plot it.
    '''
    density = gaussian_kde(exprs)
    xs = np.linspace(min(exprs), max(exprs), 200)
    density.covariance_factor = lambda : .25
    density._compute_covariance()
    plt.plot(xs, density(xs))
    plt.title('gene: %s \n mean: %s \n sd: %s' % (gene, abs_delta_mean, samp_sd))
    plt.show()

def calculate_sd_of_sample(null_mean, null_sd, sample_exprs):
    '''
    Calculate how many standard deviations away the sample_exprs
    is from the null_mean.
    
    It is an absolute number, agnostic of direction.
    '''
    '''
    print null_mean
    print null_sd
    print sample_exprs
    raw_input()
    '''
    abs_delta_mean = abs(null_mean - sample_exprs)
    if null_sd != 0:
        samp_sd = float(abs_delta_mean) / null_sd
    else:
        samp_sd = None
    return samp_sd

def get_mean_sd_from_list(null_exprs):
    '''
    From list of values, return the mean and sd of the values.
    '''
    null_exprs = [float(i) for i in null_exprs]
    # print('Null exprs: %s' %null_exprs)
    mean = float(sum(null_exprs)) / len(null_exprs)
    sd = float(np.std(null_exprs))
    return mean, sd

def both_samples_unexpressed(exprs1, exprs2, min_exprs_lvl):
    '''
    
    '''
    if exprs1 < min_exprs_lvl and exprs2 < min_exprs_lvl:
        return False
    else:
        return False

def is_outlier_expression(null_exprs, sample_exprs, sd_threshold=2.7):
    '''
    Fit null_exprs to normal distribution, extracting
    mean and sd.
    
    Compare with sample_exprs, if greater than sd_threshold, then say
    it is outlier. (True or False)
    '''
    # define constants
    min_exprs_lvl = 3
    
    null_mean, null_sd = get_mean_sd_from_list(null_exprs)
    samp_sd = calculate_sd_of_sample(null_mean, null_sd, sample_exprs)
    if samp_sd < sd_threshold or \
        both_samples_unexpressed(null_mean, sample_exprs, 
                               min_exprs_lvl=min_exprs_lvl) or \
        samp_sd == None:
        return False
    else:
        return True
    
def write_sample_gene_to_file(gene, sample, writeobj):
    '''
    Write gene, sample to writeobj, used when we know
    this gene is outlier in this sample.
    '''
    writeobj.writerow([sample, gene])
    return None

def get_exprs_from_gene_dic(null_gene_exprs_dic, gene):
    '''
    Because the function "load_gene_exprs" loads a tuple
    containing (sample name, exprs), we want to extract
    the expressions only. Removing sample names.
    '''
    samp_exprs_tup_list = null_gene_exprs_dic[gene]
    exprs_list = [i[1] for i in samp_exprs_tup_list]
    return exprs_list

def find_outlier_genes(null_gene_exprs_dic, test_samples,
                       readheader, readobj, writeobj):
    # Define column name constants
    gene_str = 'gene'
    
    # Iterate rows in readobj, assumes header already read.
    outlier_count = 0
    for row in readobj:
        gene = row[readheader.index(gene_str)]
        null_exprs = get_exprs_from_gene_dic(null_gene_exprs_dic, gene)
        # Uncomment plot_exprs to see gene exprs distribution across samp
        # plot_exprs(null_exprs, gene, mean, sd)
        for samp in test_samples:
            samp_exprs = float(row[readheader.index(samp)])
            if is_outlier_expression(null_exprs, samp_exprs):
                write_sample_gene_to_file(gene, samp, writeobj)
                outlier_count += 1
            else:
                pass
    if outlier_count % 1000 == 0:
        print('outlier_count: %s' % outlier_count)
    return outlier_count

def get_gene_exprs_distribution(null_gene_exprs_dic):
    means = []
    for gene in null_gene_exprs_dic:
        tup_list = null_gene_exprs_dic[gene]
        exprs_list = [float(i[1]) for i in tup_list]
        means.append(float(sum(exprs_list)) / len(exprs_list))
    mymean = float(sum(means)) / len(means)
    sd = np.std(means)
    plot_exprs(means, 'all', mymean, sd)

def main(exprs_file, null_samples_file, test_samples_file, output_file):
    # Define constants: output column names.
    sample_colname = 'sample'
    gene_colname = 'gene'
    
    # Get samples as list.
    null_samples = file_utilities.get_samples_from_file(null_samples_file)
    test_samples = file_utilities.get_samples_from_file(test_samples_file)
    # Create dictionary of genes:[gene_expressions]
    null_gene_exprs_dic = load_gene_exprs(exprs_file, null_samples)
    
    # To view global gene exprs distribution, uncomment get_gene_exprs_distribution
    # get_gene_exprs_distribution(null_gene_exprs_dic)
    
    # Create output objective
    jwriter = file_utilities.csv_obj(output_file, 'write')
    jreader = file_utilities.csv_obj(exprs_file, 'read')
    with jwriter:
        # Write header
        jwriter.writeobj.writerow([sample_colname, gene_colname])
        with jreader:
            readheader = jreader.readobj.next()
            outlier_count = find_outlier_genes(null_gene_exprs_dic,
                                               test_samples,
                                               readheader,
                                               jreader.readobj,
                                               jwriter.writeobj)
    print('%s rows written to file: %s' % (outlier_count, jwriter.filepath))
    
if __name__ == '__main__':
    if len(sys.argv) < 5:
        print('Gene expression file, PCa samples, NEPC samples '\
              'and output filepath must be specified in command line.')
        sys.exit()
    exprs_file = sys.argv[1]
    null_samples_file = sys.argv[2]  # e.g. PCa
    test_samples_file = sys.argv[3]  # e.g. NEPC
    output_file = sys.argv[4]
    main(exprs_file, null_samples_file, test_samples_file, output_file)
