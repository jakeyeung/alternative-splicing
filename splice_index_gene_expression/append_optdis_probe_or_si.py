'''
Created on 2013-06-04

@author: jyeung

Take a genelist from optdis output and append to each row:
1) probe information
2) probe expression/si values in each sample
'''


import sys
import csv
import os
from scipy import stats
from utilities import set_directories, read_write_gene_info

# Set constants
inputdir = 'input'
outputdir = 'output'
genesymbol_colname = 'gene_symbol'
probe_colname = 'probe'
probe_or_gene_colname = 'probe_or_gene'
sample_list = ['X946_Urethra', 'X972.2.Penila', 'AB352', 'X963.Lpa.JP', 
               'X963.L.LN', 'X1005', 'X890L', 'X945', 'X961']
genelist_colnames = [genesymbol_colname, 'subnetwork']
# last sample assumed to be last column. 


# Set directories
mydirs = set_directories.mydirs(inputdir, outputdir)

def read_probe_si_file(data_fullpath, scale_values=False, gene_exprs_table=False):
    '''
    Take probe expression or si value file and read all of its contents
    storing it into a dictionary which can easily be searchable later. 
    
    Inputs: 
    data_full_path:
    full path to either probe expression or si value file containing: 
        1) gene symbol column name
        2) probe column name
        3) column saying whether alternative splicing occured (probe or gene)
        4) first sample column to the last column are the expression or SI values.
            The expression/SI values are either standardized or not. 
    scale_values:
    Whether to scale expression values by mean centering and standard deviation normalization.
    
    gene_exprs_table:
    Whether the data_fullpath is a gene_expression_table, which is handled differently than others.
    
    Outputs: a dictionary with:
        1) Gene names as keys
        2) Gene name, probename, whether it was probe or gene, 
           and expression/SI values as values.
    '''
    errorcount = 0
    probe_si_dic = {}
    with open(data_fullpath, 'rb') as probe_si_file:
        probesi_reader = csv.reader(probe_si_file, delimiter='\t')
        probesi_headers = probesi_reader.next()
        for row in probesi_reader:
            genename = row[probesi_headers.index(genesymbol_colname)]
            if gene_exprs_table == False:
                probename = row[probesi_headers.index(probe_colname)]
                probe_or_gene = row[probesi_headers.index(probe_or_gene_colname)]
            elif gene_exprs_table == True:
                '''
                No probename or probe_or_gene columns in gene exprs table. 
                '''
                pass
            else:
                sys.exit('gene_exprs_table must be True or False')
            try:
                values = []
                for s in sample_list:
                    values.append(float(row[probesi_headers.index(s)]))
                # values = [float(i) for i in row[probesi_headers.index(first_sample):]]
                if scale_values == True:
                    '''
                    Calculate z-score of values, because we used z-scores of these 
                    values for calculating OptDis, if we are to column-bind these 
                    values, we need to make sure we're using same values as the ones
                    used to calculate subnetworks in optdis. 
                    '''
                    values = list(stats.zscore(values))
                if gene_exprs_table == False:
                    probe_si_dic[genename] = [genename, probename, 
                                              probe_or_gene] + values
                elif gene_exprs_table == True:
                    '''
                    No probename or probe_or_gene columns in gene exprs table. 
                    '''
                    probe_si_dic[genename] = [genename] + values
                else:
                    sys.exit('gene_exprs_table must be True or False')
            except ValueError:
                errorcount += 1
    # Return dictionary and all headers except the first element and 
    # unspecified samples because first element is 
    # an ensembl gene id and we have no use for it. 
    if gene_exprs_table == False:
        headers = [genename, probename, probe_or_gene] + sample_list
    elif gene_exprs_table == True:
        headers = [genename] + sample_list
    else:
        sys.exit('gene_exprs_table neither True or False, check your input.')
    print('Done, skipped %s rows probably containing NAs' %errorcount)
    return probe_si_dic, headers

def columnbind_optdis_results(expression_info_dic, read_write_object, headers):
    '''
    Purpose: take dictionary containing gene names and expression information
            columnbind the information to the read object, and save that file
            as the write object. 
    Inputs:
        expression_info_dic: genenames as keys, expression info as values. 
        read_write_object: an object for reading data and writing data, initialized 
                            by the read_write_gene_info class. 
                            The readfile contains no column 
        headers: a list containing column names to add to output.
    Outputs:
        Output written to file. 
    '''
    with read_write_object:
        # Write columns to file
        read_write_object.writenext(headers)
        # For each row, find its corresponding expression values in the
        # dictionary, then write row values AND expression values to output.
        while True:
            try:
                row = read_write_object.readnext()
            except StopIteration:
                print('%s rows read, %s rows written, breaking...'\
                       %(read_write_object.readrowcount, 
                         read_write_object.writerowcount))
                break
            gene = row[headers.index(genesymbol_colname)]
            # Find values in dictionary for this gene
            try:
                values = expression_info_dic[gene]
                read_write_object.writenext(row + values)
            except KeyError:
                print('%s not found, skipping...' %gene)
    
            
if __name__ == '__main__':
    if len(sys.argv) < 5:
        print('Probe/SI info path, true/false for gene_exprs_table, optdis results path '\
              ' and output filename must be provided in command line.')
        sys.exit()
    '''
    Example inputs:
    sys.argv[1] = tables/si_and_gene_expression_data.txt
    sys.argv[2] = optdis_outputs/si_and_gene_exprs/Results/MarkerDiscovery/
                    takeda/subtype/INFORGAIN/full_list_si_and_gene_exprs.txt
    sys.argv[3] = optdis_outputs/si_and_gene_exprs/Results/MarkerDiscovery/
                    takeda/subtype/INFORGAIN/full_list_si_and_gene_exprs_appended.txt
    '''
    data_partialpath = sys.argv[1]
    gene_exprs_table = sys.argv[2]
    optdis_partialpath = sys.argv[3]
    output_partialpath = sys.argv[4]
    
    if gene_exprs_table in ['true', '1', 't', 'y', 'yes', 'True', 'TRUE']:
        gene_exprs_table = True
    elif gene_exprs_table in ['False', 'FALSE', 'F', 'f', 'no', 'No', 'NO']:
        gene_exprs_table = False
    else:
        sys.exit('Unknown input %s, must be either True or False.' %gene_exprs_table)

    data_fullpath = os.path.join(mydirs.outputdir, data_partialpath)
    optdis_fullpath = os.path.join(mydirs.outputdir, optdis_partialpath)
    output_fullpath = os.path.join(mydirs.outputdir, output_partialpath)
    
    probe_si_dic, headers = read_probe_si_file(data_fullpath, 
                                               scale_values=True, 
                                               gene_exprs_table=gene_exprs_table)
    
    read_write_optdisresults = read_write_gene_info.read_write_gene_info(optdis_fullpath, 
                                                                         output_fullpath, 
                                                                         header=False)
    
    # Add gene list headers
    headers = genelist_colnames + headers
    columnbind_optdis_results(probe_si_dic, read_write_optdisresults, headers)
    
    
    
    
    
    
    