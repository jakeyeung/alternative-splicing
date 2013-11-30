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

def get_output_colnames(collapsed=False):
    '''
    Get column names... simple as that.
    
    Collapsed version is reduced (no pattern name)
    '''
    # def constants
    pattern_name_str = 'pattern_name'
    gene_name_str = 'rbp_name'
    region_str = 'exon_intron_region'
    occurences_str = 'n_motif_occurences'
    med_qval_str = 'median_q_value'
    incl_excl_str = 'inclusion_or_exclusion'
    if collapsed == False:
        return [pattern_name_str, gene_name_str, region_str, 
                occurences_str, med_qval_str, incl_excl_str]
    elif collapsed == True:
        return [gene_name_str, region_str, 
                occurences_str, med_qval_str, incl_excl_str]
    else:
        print 'Collapsed should be True or False. %s found.' %collapsed
        sys.exit()
    
def get_fimo_subkeys():
    '''
    Get fimo subkeys to store into subdictionary.
    Subkeys:
    rbp_name
    exon_intron_region
    inclusion_or_exclusion
    
    Does not include q-value...
    Return as a list.
    '''
    rbpname = 'rbp_name'
    region = 'exon_intron_region'
    incl_excl = 'inclusion_or_exclusion'
    return [rbpname, region, incl_excl]
    
def get_info_from_fimo_output(fimo_path, writeobj, fimo_dic, 
                              remove_repeats=False):
    '''
    Read fimo output, it's not ordered in any sane way, so we will do it 
    with no shortcuts.
    
    Create a dictionary, keys will be pattern name,
    values will be list of qvalues,
    
    Store into dictionary the following:
        pattern name (easier to retrieve later)
        median qvalue
        n_occurences
    
    Remove repeats: True by default. It will not count rows
    of same pattern name and same sequence name as previous row.
    '''
    # Def colname constants
    pattern_name_colname = '#pattern name'
    qval_colname = 'q-value'
    seq_name_colname = 'sequence name'
    pattern_name_subkey = 'pattern_name'
    qval_med_subkey = 'median_q_value'
    n_occurs_subkey = 'n_motif_occurences'
    
    with open(fimo_path, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        header = myreader.next()
        for row in myreader:
            '''
            Check every row to make sure seq name and pattern name
            are different. If both are the same, then we don't want to
            include it in our analysis (too many RBPs)
            '''
            seq_name = row[header.index(seq_name_colname)]
            pattern_name = row[header.index(pattern_name_colname)]
            qval = row[header.index(qval_colname)]
            # Create new key for pattern name if it doesnt exist in dic
            # Create qval subkey with empty list.
            if pattern_name not in fimo_dic:
                fimo_dic[pattern_name] = {}
                fimo_dic[pattern_name][qval_colname] = []
                fimo_dic[pattern_name][seq_name_colname] = []
            # Check that seq name does not already exist in our dic.
            # if it is already in our dic, then move to next row.
            if seq_name in fimo_dic[pattern_name][seq_name_colname]:
                continue
            # Append empty list with qvalues from each row...
            fimo_dic[pattern_name][qval_colname].append(qval)
            fimo_dic[pattern_name][seq_name_colname].append(seq_name)
    # Get median q-value, n_motif_occurences and pattern name to dictionary.
    for pat_name, qval_dic in fimo_dic.iteritems():
        # Get pattern name into dictionary info
        fimo_dic[pat_name][pattern_name_subkey] = pat_name
        # calculate med qval
        qval_list = [float(i) for i in qval_dic[qval_colname]]
        med_qval = stats_functions.median(qval_list)
        # Store med qval to dic
        fimo_dic[pat_name][qval_med_subkey] = med_qval
        # Get n_occurences
        fimo_dic[pat_name][n_occurs_subkey] = len(qval_list)
    return fimo_dic

def get_number_of_sequences_from_fimohtml(fimo_html_path, rownumber=42):
    '''
    The 42nd line from fimo html output contains total number
    of fasta sequences used to run FIMO.
    '''
    currentrow = 0
    with open(fimo_html_path, 'rb') as readfile:
        while currentrow != rownumber:
            rowstr = readfile.readline()
            currentrow += 1
        # Strip of any white space
        total_sequences = float(rowstr.strip())
    return total_sequences

def get_info_from_fimo_output2(fimo_path, writeobj, fimo_dic, 
                               convert_to_fraction=False):
    '''
    Second attempt at getting info from fimo.
    I want gene name to show up only ONCE for a given sequence.
    '''
    # Def colname constants
    pattern_name_colname = '#pattern name'
    qval_colname = 'q-value'
    seq_name_colname = 'sequence name'
    pattern_name_subkey = 'pattern_name'
    qval_med_subkey = 'median_q_value'
    n_occurs_subkey = 'n_motif_occurences'
    
    with open(fimo_path, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        header = myreader.next()
        for row in myreader:
            '''
            Check every row to make sure seq name and pattern name
            are different. If both are the same, then we don't want to
            include it in our analysis (too many RBPs)
            '''
            seq_name = row[header.index(seq_name_colname)]
            pattern_name = row[header.index(pattern_name_colname)]
            qval = row[header.index(qval_colname)]
            '''
            # Get gene name from pattern name, expect 
            string to be: GENENAME,MOTIF,IorD, so we will separate commas.
            '''
            gene_name = pattern_name.split(',')[0]
            
            # Create new key for pattern name if it doesnt exist in dic
            # Create qval subkey with empty list.
            if gene_name not in fimo_dic:
                fimo_dic[gene_name] = {}
                fimo_dic[gene_name][qval_colname] = []
                fimo_dic[gene_name][seq_name_colname] = []
            # Check that seq name does not already exist in our dic.
            # if it is already in our dic, then move to next row.
            if seq_name in fimo_dic[gene_name][seq_name_colname]:
                continue
            # Append empty list with qvalues from each row...
            fimo_dic[gene_name][qval_colname].append(qval)
            fimo_dic[gene_name][seq_name_colname].append(seq_name)
    # Get median q-value, n_motif_occurences and pattern name to dictionary.
    for gene_name, qval_dic in fimo_dic.iteritems():
        # Get pattern name into dictionary info
        fimo_dic[gene_name][pattern_name_subkey] = gene_name
        # calculate med qval
        qval_list = [float(i) for i in qval_dic[qval_colname]]
        med_qval = stats_functions.median(qval_list)
        # Store med qval to dic
        fimo_dic[gene_name][qval_med_subkey] = med_qval
        # Get n_occurences
        # fimo_dic[gene_name][n_occurs_subkey] = len(qval_list)
        n_occurences = len(qval_list)
        if convert_to_fraction == True:
            '''
            # Get n_occurences divided by total fasta inputs
            # Get dirname of fimo txt file, append the html file
            # to read html file.
            '''
            fimo_html_path = os.path.join(os.path.dirname(fimo_path), 'fimo.html')
            n_fasta_seqs = \
                get_number_of_sequences_from_fimohtml(fimo_html_path, rownumber=42)
            fimo_dic[gene_name][n_occurs_subkey] = float(n_occurences) / n_fasta_seqs
        elif convert_to_fraction == False:
            '''
            Just get the occurences, no dividing by fasta seqs
            '''
            fimo_dic[gene_name][n_occurs_subkey] = n_occurences
        else:
            print 'Expected convert to fraction '\
            'to be True or False. %s found.' %convert_to_fraction
    return fimo_dic

def get_rbp_name_from_pat_name(pat_name):
    '''
    From a pat name of expected CSV form: rbp,motifid,D_or_I
    extract rbp.
    '''
    pat_name_split = pat_name.split(',')
    # Get first element, rejoin...
    rbp = pat_name_split[0]
    return rbp

def add_info_for_dic(fimo_info_dic, region, incl_or_excl):
    '''
    Inputs:
    fimo_info_dic: dic containing keys of rbp pattern name and med qval.
    region: exon_intron region, like exon_1 or intron_2_3p
    incl_or_excl: inclusion or exclusion.
    Add rbp_name, region and incl_or_excl to the dictionary.
    '''
    for pat_name in fimo_info_dic.keys():
        rbp_name = get_rbp_name_from_pat_name(pat_name)
        subkeys = get_fimo_subkeys()
        for subk, subval in zip(subkeys, [rbp_name, region, incl_or_excl]):
            fimo_info_dic[pat_name][subk] = subval
    return fimo_info_dic

def add_info_for_dic2(fimo_info_dic, region, incl_or_excl):
    '''
    Second try: now pat name is a gene name. No need to get rbp_name from
    a function.
    '''
    for rbp_name in fimo_info_dic.keys():
        subkeys = get_fimo_subkeys()
        for subk, subval in zip(subkeys, [rbp_name, region, incl_or_excl]):
            fimo_info_dic[rbp_name][subk] = subval
    return fimo_info_dic

def get_region_from_dirname(mydir):
    '''
    Parse name of the directory and try to get
    the region (exon_1, incl_2_3p...)
    
    Expecting mydir to be of form
    ExonOrIntron_1or2or3_InclusionOrExclusion
    '''
    mydir_split = mydir.split('_')
    # Remove last element, keep the first two elements, rejoin.
    mydir = '_'.join(mydir_split[:-1])
    return mydir

def get_incl_or_excl_from_dirname(mydir):
    '''
    Parse name of the directory and try to get
    whether it is inclusion or exclusion...
    '''
    mydir_split = mydir.split('_')
    # keep last element only, rejoin.
    incl_or_excl = mydir_split[-1]
    if incl_or_excl in ['inclusion', 'exclusion', 
                        'inclusion.shuffled', 'exclusion.shuffled']:
        return incl_or_excl
    else:
        print 'Expected incl_or_excl to be "inclusion" or "exclusion". '\
        '%s found.' %incl_or_excl
        sys.exit()
        
def write_fimo_dic_to_file(mydic, mywriter, collapsed):
    '''
    After reading FIMO output into a dictionary for
    a single region (e.g. exon_1, let's write what we have
    to the output file.)
    We want column order to be (from left to right) the same as 
    get_output_colnames(), see that function for current order of
    columns.
    
    As of 20-11-2013, the colnames are:
    pattern_name_str = 'pattern_name'
    gene_name_str = 'rbp_name'
    region_str = 'exon_intron_region'
    occurences_str = 'n_motif_occurences'
    med_qval_str = 'median_q_val'
    incl_excl_str = 'inclusion_or_exclusion'
    '''
    colnames = get_output_colnames(collapsed=collapsed)
    writecount = 0
    for pattern_name in mydic.keys():
        writerow = []    # init
        # Get information according to colnames. 
        # Colnames should match subkeys in dic
        for subkey in colnames:
            # Assumes values in subkey are not lists...
            writerow.append(mydic[pattern_name][subkey])
        mywriter.writerow(writerow)
        writecount += 1
    return writecount

def reduce_dic_rbps_only(fimo_info_dic):
    '''
    Reduce dic to have ONE rbp per dic.
    n_motif_occurences will then be the sum of all the rbp_names.
    '''
    reduced_dic = {}
    n_motif_subkey = 'n_motif_occurences'
    incl_excl_subkey = 'inclusion_or_exclusion'
    exon_intron_subkey = 'exon_intron_region'
    med_q_subkey = 'median_q_value'
    rbp_name_subkey = 'rbp_name'
    for pat_name in fimo_info_dic.keys():
        rbp_name = fimo_info_dic[pat_name][rbp_name_subkey]
        if rbp_name not in reduced_dic:
            '''
            # Init values for n motif subkey
            # incl/excl, exon/intron, do not change, so define them once 
            (by defining them within this if)
            Median q values are not too useful, just put them as med_q_subkey, 
            but I wouldnt take them to heart.
            '''
            reduced_dic[rbp_name] = {}
            reduced_dic[rbp_name][n_motif_subkey] = 0
            for k in [incl_excl_subkey, exon_intron_subkey, 
                      med_q_subkey, rbp_name_subkey]:
                reduced_dic[rbp_name][k] = fimo_info_dic[pat_name][k]
        reduced_dic[rbp_name][n_motif_subkey] += \
            fimo_info_dic[pat_name][n_motif_subkey]
    return reduced_dic

def convert_to_boolean(myvariable):
    '''
    Expects True or False as input as a string.
    Outputs True or False as boolean variable.
    '''
    if myvariable in ['True', 'TRUE', 'true', 't', 'T']:
        myvariable = True
    elif myvariable in ['False', 'FALSE', 'false', 'f' , 'F']:
        myvariable = False
    else:
        print 'Expected variable to be True or False. %s found.' \
            %myvariable
        sys.exit() 
    return myvariable

def main():
    usage = 'usage: %prog fimo_out_dir output_file.txtfile\n'\
        'Two args must be specified in commandline: \n'\
        '1) directory containing exon/intron-region directories, '\
        'inside which contains fimo outputfiles.\n'\
        '2) Output file.\n'
    parser = OptionParser(usage=usage)
    parser.add_option('-f', '--fimo_filename', dest='fimo_filename',
                      help='Name of fimo output filename. Default fimo.txt',
                      default='fimo.txt')
    parser.add_option('-c', '--collapse_rbps', dest='collapse_rbps',
                      help='Collapse all motifs for one rbp into one line. Default True',
                      default='True')
    parser.add_option('-C', '--convert_to_fraction', dest='convert_to_fraction',
                      help='Convert number of hits to fraction of total sequences. Default True',
                      default='True')
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print(usage)
        sys.exit()
    fimo_dir = args[0]
    out_path = args[1]
    fimo_filename = options.fimo_filename
    
    # convert True/False options to Boolean
    collapse_rbps = convert_to_boolean(options.collapse_rbps)
    convert_to_fraction = convert_to_boolean(options.convert_to_fraction)
    '''
    if options.collapse_rbps in ['True', 'TRUE', 'true', 't', 'T']:
        collapse_rbps = True
    elif options.collapse_rbps in ['False', 'FALSE', 'false', 'f' , 'F']:
        collapse_rbps = False
    else:
        print 'Expected collapse_rbps to be True or False. %s found.' \
            %collapse_rbps
        sys.exit() 
    '''
    
    # Check that outpath does not already exist. If it does, then abort.
    if os.path.isfile(out_path):
        print '%s already exists. Aborting FIMO output parsing.' \
            %out_path
        sys.exit() 
    
    # Get list of exon-intron directories (exon1, exon2... intron2_3p)
    # I only want directories, so check that it is a directory, not a file.
    dirs = \
        [d for d in os.listdir(fimo_dir) if \
            os.path.isdir(os.path.join(fimo_dir, d))]
    
    # Initialize writefile
    writecount = 0
    with open(out_path, 'wb') as writefile:
        mywriter = csv.writer(writefile, delimiter='\t')
        # Write header
        header = get_output_colnames(collapsed=collapse_rbps)
        mywriter.writerow(header)
        # Loop through directories
        for d in dirs:
            fimo_info_dic = {}
            region = get_region_from_dirname(d)
            incl_or_excl = get_incl_or_excl_from_dirname(d)
            
            fimo_path = os.path.join(fimo_dir, d, fimo_filename)
            if collapse_rbps:
                fimo_info_dic = \
                    get_info_from_fimo_output2(fimo_path, 
                                               mywriter, 
                                               fimo_info_dic,
                                               convert_to_fraction=\
                                                    convert_to_fraction)
                fimo_info_dic = add_info_for_dic2(fimo_info_dic, region, 
                                                  incl_or_excl)
            else:
                fimo_info_dic = get_info_from_fimo_output(fimo_path, 
                                                          mywriter, 
                                                          fimo_info_dic)
                fimo_info_dic = add_info_for_dic(fimo_info_dic, region, 
                                                 incl_or_excl)
            '''
            Loop through keys and do a double check that there are no sequence names duplicated.
            '''
            for k in fimo_info_dic.keys():
                if len(fimo_info_dic[k]['sequence name']) != \
                    len(list(set(fimo_info_dic[k]['sequence name']))):
                    print 'Warning: %s has duplicate sequence names. '\
                    'Please double check.' %k
                    raw_input('Press enter to continue...')
            '''
            Reduce dic to include only RBP (collapse different motifs 
            from same RBP into one)
            '''
            if collapse_rbps:
                fimo_info_dic = reduce_dic_rbps_only(fimo_info_dic)
            # Write dic to file..
            writecount += write_fimo_dic_to_file(fimo_info_dic, mywriter, 
                                                 collapsed=collapse_rbps)
    print '%s rows written to file: %s' %(writecount, out_path)

if __name__ == '__main__':
    main()