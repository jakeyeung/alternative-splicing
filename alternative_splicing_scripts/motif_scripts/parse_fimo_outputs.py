'''
Created on 2013-11-15

@author: jyeung
FIMO a textfile at each exon/intron region.
Parse them all, and put them in a textfile that can be
easily ggplot'd.
'''

import sys
import os
import csv
from optparse import OptionParser
from utilities import stats_functions

def get_output_colnames():
    '''
    Get column names... simple as that.
    '''
    # def constants
    pattern_name_str = 'pattern_name'
    gene_name_str = 'gene_name'
    region_str = 'exon_intron_region'
    occurences_str = 'n_motif_occurences'
    med_qval_str = 'median_q_val'
    incl_excl_str = 'inclusion_or_exclusion'
    return [pattern_name_str, gene_name_str, region_str, 
            occurences_str, med_qval_str, incl_excl_str]
    
def get_fimo_subkeys():
    '''
    Get fimo subkeys to store into subdictionary.
    Subkeys:
    q-value
    rbp_name
    n_motif_occurences
    exon_intron_region
    inclusion_or_exclusion
    
    Return as a list.
    '''
    qval = 'q-value'
    rbpname = 'rbp_name'
    n_occurs = 'n_motif_occurences'
    region = 'exon_intron_region'
    incl_excl = 'inclusion_or_exclusion'
    return [qval, rbpname, n_occurs, region, incl_excl]
    
def get_qvals_from_fimo(fimo_path, writeobj):
    '''
    Read fimo output, it's not ordered in any sane way, so we will do it 
    with no shortcuts.
    
    Create a dictionary, keys will be pattern name,
    values will be list of qvalues,
    '''
    # Def colname constants
    pattern_name_colname = '#pattern name'
    qval_colname = 'q-value'
    qval_med_colname = 'median_q-value'
    
    fimo_dic = {}
    with open(fimo_path, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        header = myreader.next()
        for row in myreader:
            pattern_name = row[header.index(pattern_name_colname)]
            qval = row[header.index(qval_colname)]
            # Create new key for pattern name if it doesnt exist in dic
            # Create qval subkey with empty list.
            if pattern_name not in fimo_dic:
                fimo_dic[pattern_name] = {}
                fimo_dic[pattern_name][qval_colname] = []
            # Append empty list with qvalues from each row...
            fimo_dic[pattern_name][qval_colname].append(qval)
    # Get median q-value for each list...
    for pat_name, qval_list in fimo_dic.iteritems():
        # Get median qvalue
        med_qval = stats_functions.median(qval_list)
        fimo_dic[pat_name][qval_med_colname] = med_qval
    return fimo_dic

def add_info_for_dic(fimo_info_dic, region, incl_or_excl):
    '''
    Inputs:
    fimo_info_dic: dic containing keys of rbp pattern name and med qval.
    region: exon_intron region, like exon_1 or intron_2_3p
    incl_or_excl: inclusion or exclusion.
    '''
    pass

def get_region_from_dirname(mydir):
    '''
    Parse name of the directory and try to get
    the region (exon_1, incl_2_3p...)
    '''
    pass

def get_incl_or_excl_from_dirname(mydir):
    '''
    Parse name of the directory and try to get
    whether it is inclusion or exclusion...
    '''
    pass

def main():
    usage = 'usage: %prog fimo_out_dir output_file.txtfile\n'\
        'Two args must be specified in commandline: \n'\
        'directory containing exon/intron-region directories, '\
        'inside which contains fimo outputfiles.\n'\
        'Output file.\n'
    parser = OptionParser(usage=usage)
    parser.add_option('-f', '--fimo_filename', dest='fimo_filename',
                      help='Name of fimo output filename. Default fimo.txt',
                      default='fimo.txt')
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print(usage)
        sys.exit()
    fimo_dir = args[0]
    out_path = args[1]
    fimo_filename = options.fimo_filename
    
    # Get list of exon-intron directories (exon1, exon2... intron2_3p)
    dirs = os.listdir(fimo_dir)
    
    # Initialize writefile
    with open(out_path, 'wb') as writefile:
        mywriter = csv.writer(writefile, delimiter='\t')
        # Write header
        header = get_output_colnames()
        mywriter.writerow(header)
        # Loop through directories
        for d in dirs:
            region = get_region_from_dirname(d)
            incl_or_excl = get_incl_or_excl_from_dirname(d)
            
            fimo_path = os.path.join(fimo_dir, d, fimo_filename)
            fimo_info_dic = get_qvals_from_fimo(fimo_path, mywriter)
            fimo_info_dic = add_info_for_dic(fimo_info_dic, region, 
                                             incl_or_excl)
        

if __name__ == '__main__':
    main()