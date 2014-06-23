'''
Created on 2014-02-21

@author: jyeung

Run anchor in batch mode.

Reads as input the output from create_dna_protein_summary_file.py

Reads file and creates a fasta file for each sequence, uses that
created fasta file as input to ANCHOR program.

Does this in parallel fashion, creating a large number of files
and then reading those files and summarizing the data.

Any more than 500 files, and the program should complain
unless explicitly forced to run.

Must run on a machine with ANCHOR installed, with PATH and
ANCHOR_PATH variables set appropriately.

e.g. ANCHOR_PATH=/home/jyeung/ANCHOR
    PATH=$PATH:/home/jyeung/ANCHOR
'''

import sys
import os
import csv
from optparse import OptionParser
from utilities.anchor import AnchorOutput

class Printer():
    """
    Print things to stdout on one line dynamically
    """

    def __init__(self,data):

        sys.stdout.write("\r\x1b[K"+data.__str__())
        sys.stdout.flush()

def get_dic_from_protein_file(protein_file,
                              gene_colname,
                              event_colname,
                              rf_colname,
                              seq_colname):
    '''
    Read protein summary file and index information into
    a dictionary.

    Dictionary is of form:
    {event_rf: {gene_colname: gene, seq_colname: seq}}

    event_rf is concatenation between gene and event.
    '''
    protein_dic = {}
    with open(protein_file, 'rb') as readfile:
        myreader = csv.reader(readfile, delimiter='\t')
        header = myreader.next()
        for row in myreader:
            gene = row[header.index(gene_colname)]
            event = row[header.index(event_colname)]
            # replace : with _ in events, otherwise it
            # doesnt work well with windows filenames
            event = event.replace(':', '_')
            rf = row[header.index(rf_colname)]
            seq = row[header.index(seq_colname)]
            # create dic key from event and rf
            key = '_'.join([event, rf])
            if key not in protein_dic:
                protein_dic[key] = {}
                for subkey, subval in \
                    zip([gene_colname, seq_colname], [gene, seq]):
                    protein_dic[key][subkey] = subval
            else:
                print 'Unexpected duplicate key: %s' %key
                sys.exit()
    return protein_dic

def create_anchor_input_file(anchor_dir, event_fname, seq, ext):
    '''
    Creates an anchor input file from input sequence.
    The input filename is event_fname.ext

    Inputs:
        anchor_dir: directory where file will be created.
        event_fname: usually event_id:reading_frame
        seq: amino acid sequence to be used
        ext: filename extension, e.g. anchorinput

    Output:
        anchor input filename called:
            event_fname.ext
            in directory:
                anchor_dir
        in which amino acid sequence is inside.
    '''
    file_fullname = os.path.join(anchor_dir, '.'.join([event_fname, ext]))
    with open(file_fullname, 'wb') as writefile:
        # write two lines: header and sequence, end with new line.
        writefile.write(''.join(['>', event_fname, '\n']))
        writefile.write(''.join([seq, '\n']))
    return file_fullname

def run_anchor(input_file, output_file):
    '''
    Opens a bash subprocess and runs ANCHOR.

    Anchor command is: anchor inputfile > outputfile

    Needs PATH and ANCHOR_PATH variables set already.
    '''
    anchor_command = ' '.join(['anchor', input_file, '>', output_file])
    os.system(anchor_command)
    # check output file was created
    if not os.path.exists(output_file):
        print '%s did not get created.' %output_file
        #print 'Error in running command:\n%s\nSkipping...' %anchor_command
        return None
    return output_file

def get_basename(mypath, remove_ext=True):
    '''
    Given path, returns base name. Removing extension is an option.
    '''
    base_with_ext = os.path.basename(mypath)
    if remove_ext:
        base = os.path.splitext(base_with_ext)[0]
        return base
    else:
        return base_with_ext

def get_corresponding_output(input_path, output_dir, ext):
    '''
    Given input path, retrieve basename, replace with
    user-defined extension.

    Then join output_dir with new extension'd filename.
    '''
    base = get_basename(input_path, remove_ext=True)
    # add new extension
    base_out_ext = '.'.join([base, ext])
    # join new path
    output_path = os.path.join(output_dir, base_out_ext)
    return output_path

def write_dic_to_file(outdic, summary_file,
                      event_colname,
                      gene_colname,
                      rf_colname,
                      seq_colname,
                      binding_regions_colname):
    '''
    Write outdic to summary file.

    Outdic format:
    {ID: {amino_acid_sequence: MySeq, gene_name:
    MyGene: binding:regions:{subdic}}}

    Where ID is like:
    chr19_5229502_5229695_-@chr19_5229327_5229353_-@chr19_5225738_5225855_-_0

    Replace '_' with ':' in ID, and split off the last integer as reading frame
    to recreate event_id and reading frame separately.

    Write to file with column names:
    event_colname, gene_colname, rf_colname, seq_colname, binding_regions
    '''
    colnames = [event_colname, gene_colname,
                rf_colname, seq_colname,
                binding_regions_colname]

    with open(summary_file, 'wb') as writefile:
        mywriter = csv.writer(writefile, delimiter='\t')
        mywriter.writerow(colnames)    # header
        for rowcount, custom_id in enumerate(outdic):
            # retrieve event_id and reading_frame from custom_id
            event_id = custom_id.split('_')[:-1]    # last element is RF
            event_id = ':'.join(event_id)    # recreate original ID
            rf = custom_id.split('_')[-1]
            gene = outdic[custom_id][gene_colname]
            seq = outdic[custom_id][seq_colname]
            binding_regions = outdic[custom_id][binding_regions_colname]
            # write in same order as header
            mywriter.writerow([event_id, gene, rf, seq, binding_regions])
    print '%s rows written to: %s' %(rowcount, summary_file)
    return None

def main():
    usage = 'usage: %prog [opt] protein_summary_file output_dir output_file'\
        '\nThree arguments must be specified in command line:\n'\
        '1) protein summary file from create_dna_protein_summary_file.py\n'\
        '2) output_dir where anchor input/output files will be created. \n'\
        '3) output file, located in output directory.\n'
    parser = OptionParser(usage=usage)
    parser.add_option('--max_files', dest='max_files',
                      default=1000,
                      help='Max files to create, any larger '\
                        'and program complains')
    parser.add_option('--gene_colname', dest='gene_colname',
                      default='gene_name',
                      help='Colname containing gene name in input file')
    parser.add_option('--event_colname', dest='event_colname',
                      default='miso_event',
                      help='Colname containing event id in input file')
    parser.add_option('--rf_colname', dest='rf_colname',
                      default='reading_frame',
                      help='Reading frame colname in input file')
    parser.add_option('--seq_colname', dest='seq_colname',
                      default='amino_acid_sequence',
                      help='Amino acid sequence colname in input file')
    parser.add_option('--anchor_input_ext', dest='in_ext',
                      default='anchorinput',
                      help='Extension for anchor input files to be created.')
    parser.add_option('--anchor_output_ext', dest='out_ext',
                      default='anchoroutput',
                      help='Extension for anchor output file. '\
                        'Default anchoroutput')
    (options, args) = parser.parse_args()

    if len(args) != 3:
        print 'Requires 3 arguments to be specified in command line'
        print usage
        sys.exit()
    # parse args
    protein_file = args[0]
    output_dir = args[1]
    summary_file = args[2]
    # parse ops
    max_files = int(options.max_files)
    gene_colname = options.gene_colname
    event_colname = options.event_colname
    rf_colname = options.rf_colname
    seq_colname = options.seq_colname
    in_ext = options.in_ext
    out_ext = options.out_ext
    binding_regions_colname = 'binding_regions'

    # create output_dir if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print 'Created directory: %s' %output_dir

    # read protein file, put all sequences into a dictionary
    protein_dic = get_dic_from_protein_file(protein_file,
                                            gene_colname=gene_colname,
                                            event_colname=event_colname,
                                            rf_colname=rf_colname,
                                            seq_colname=seq_colname)
    if len(protein_dic.keys()) >= max_files:
        # complain that there are too many files to prevent
        # creating too many files.
        print 'Number of files to create (%s) exceed maximum recommended (%s). '\
            '\nIncrease --max_files to override this complaint.'
        sys.exit()

    # define anchor input directory and create it if necessary
    anchor_input_dir = os.path.join(output_dir, 'anchor_inputs')
    if not os.path.exists(anchor_input_dir):
        os.makedirs(anchor_input_dir)
        print 'Created anchor input directory: %s' %anchor_input_dir

    # define anchor output directory and create it if necessary
    anchor_output_dir = os.path.join(output_dir, 'anchor_outputs')
    if not os.path.exists(anchor_output_dir):
        os.makedirs(anchor_output_dir)
        print 'Created anchor output directory: %s' %anchor_output_dir

    print 'Creating %s files for ANCHOR...' %len(protein_dic.keys())
    # create input files for each amino acid sequence
    anchor_inputs = []
    for event_fname, subdic in protein_dic.iteritems():
        seq = subdic[seq_colname]
        anchor_inputs.append(create_anchor_input_file(anchor_input_dir,
                                                      event_fname,
                                                      seq,
                                                      ext=in_ext))
    n_inputs = len(anchor_inputs)
    print 'Created %s input files for ANCHOR.' %n_inputs

    anchor_count = 0
    failcount = 0
    anchor_outputs = []
    for input_file in anchor_inputs:
        # create corresponding outputfile, store in list for later retrieval
        output_file = get_corresponding_output(input_file,
                                               anchor_output_dir,
                                               ext=out_ext)
        output_file = run_anchor(input_file, output_file)
        if output_file is not None:    # if None: anchor failed to run.
            anchor_outputs.append(output_file)
            anchor_count += 1
            output = 'ANCHOR outputs: %s/%s' %(anchor_count, n_inputs)
            Printer(output)
        else:
            failcount += 1
    print '\n'    # adds carriage return to Printer output only after loop.
    if failcount > 0:
        # inform user of any failures if any.
        print '%s files failed to run anchor.' %failcount

    # Parse ANCHOR outputs
    for output_file in anchor_outputs:
        anchor = AnchorOutput(output_file)
        binding_regions_dic = anchor.binding_regions()
        # append to existing dic
        key = get_basename(output_file, remove_ext=True)
        protein_dic[key][binding_regions_colname] = binding_regions_dic

    # Write to output file
    summary_path = os.path.join(output_dir, summary_file)
    write_dic_to_file(protein_dic, summary_path,
                      event_colname,
                      gene_colname,
                      rf_colname,
                      seq_colname,
                      binding_regions_colname)

if __name__ == '__main__':
    main()
