'''
Created on 2013-12-19

@author: jyeung

Code adotped from cbayly
'''


import pickle
import sys 
import csv

def get_chromosome_list(add_chr_prefix=False):
    '''
    Get chromosome names, 1 2 3 4 .. 22 X Y
    
    If you want chr1, chr2...  chrX, chrY, use
    add_prefix = True, otherwise it's False
    '''
    chr_list = [str(i) for i in range(1, 23)] + ['X', 'Y']
    
    if add_chr_prefix == True:
        chr_list = [''.join(['chr', chr_name]) for chr_name in chr_list]
        
    return chr_list

def get_subkeys():
    '''
    For subkeys, we want to easily retrieve it by a function.
    
    That way we can access them without hardcoding constants.
    '''
    # def subkey names
    reading_frame_str = 'reading_frame'
    gene_name_str = 'gene_name'
    exon_number_str = 'exon_number'
    gene_id_str = 'gene_id'
    transcript_id_str = 'transcript_id'
    exon_id_str = 'exon_id'
    return reading_frame_str, gene_name_str, exon_number_str, \
            gene_id_str, transcript_id_str, exon_id_str
            
def extract_info_from_exon(exon_info, search=''):
    '''
    the 9th column in ensembl GTF file contains semi-colon separated
    values that contains gene id, gene name, etc...
    
    Input the 9th column (as a string) to this function, ask for what to search,
    and it will output the gene info you're searching for, removing any quotes
    it may be surrounded by.
    
    exon_info_as_list looks like:
    ['gene_id', '"ENSG00000142449"', 'transcript_id', '"ENST00000601739"', 'exon_num
    ber', '"46"', 'gene_name', '"FBN3"', 'gene_biotype', '"protein_coding"', 'transc
    ript_name', '"FBN3-002"', 'exon_id', '"ENSE00001663237"']
    
    Find index of search string, return the value of index + 1
    '''
    # remove semi colons, then split by space
    exon_info_as_list = exon_info.replace(';', '').split(' ')
    # search string of interest, find its index, 
    # then return the value at index + 1
    try:
        index_of_interest = exon_info_as_list.index(search)
    except ValueError:
        print 'Could not find %s in %s' %(search, exon_info)
        sys.exit()
    info_of_interest = exon_info_as_list[index_of_interest + 1]
    # Remove quotes from info of interest
    info_of_interest = info_of_interest.replace('"', '')
    return info_of_interest

def make_dictionary(ensdatabase, ensdictionary):
    ensdict = {}
    f = open(ensdatabase, 'rb')            #inputfile1 = '100part_Homo_sapiens.GRCH37.72.gtf'
    #    Grabs readframe, constructs chromosome location, constructs exon 'identifiers' (transcriptnumber, jtype(like "exon"), number).
    
    # Define index names
    chr_index = 0    # chromosome, we want 1, 2, 3.. X and Y's only
    # feature_index = 2    # we want coding sequence only
    exon_start_index = 3
    exon_end_index = 4
    strand_index = 6
    reading_frame_index = 7
    exon_info_index = 8
    
    reading_frame_str, \
    gene_name_str, \
    exon_number_str, \
    gene_id_str, \
    transcript_id_str, \
    _ = get_subkeys()
    
    chr_list = get_chromosome_list()
    
    jreader = csv.reader(f, delimiter='\t')
    
    count = 0
    for row in jreader:
        chr_name = row[chr_index]
        if chr_name not in chr_list:
            continue
        else:
            reading_frame = row[reading_frame_index]
            if reading_frame == '.':
                continue
            else:
                # Build key with chr:start:end:strand
                start = row[exon_start_index]
                end = row[exon_end_index]
                strand = row[strand_index]
                
                # Add chr to chr_name, 
                # then create location chr:start:end:strand
                # This matches MISO coordinates for cassette exons
                chr_name_prefixed = ''.join(['chr', chr_name])
                location = ':'.join([chr_name_prefixed, start, end, strand])
                
                # Get exon info from ensembl
                exon_info = row[exon_info_index]
                
                # Parse exon info to get gene name, exon numb, and ids
                gene_name = \
                    extract_info_from_exon(exon_info, search=gene_name_str)
                exon_number = \
                    extract_info_from_exon(exon_info, search=exon_number_str)
                gene_id = \
                    extract_info_from_exon(exon_info, search=gene_id_str)
                transcript_id = \
                    extract_info_from_exon(exon_info, search=transcript_id_str)
                # if location does not exist, then initialize with empty dictionary,
                # then add subdic information as a list.
                # if location exists
                if location not in ensdict:
                    # location not yet exist, we need to 
                    # initialize a subdic
                    # then initialize an empty list in each subkey
                    # before we can append as a list.
                    ensdict[location] = {}
                    # init my loops
                    for subkey, subvalue in \
                    zip([reading_frame_str, gene_name_str, exon_number_str, 
                         gene_id_str, transcript_id_str], 
                        [reading_frame, gene_name, exon_number, gene_id, 
                         transcript_id]):
                        # since new dic, add empty list, then append
                        ensdict[location][subkey] = []
                else:
                    # location exists,
                    # I assume subvalues are in list form, do nothing
                    # so I can simply 
                    # append subvalues to subkeys
                    pass
                for subkey, subvalue in \
                zip([reading_frame_str, gene_name_str, exon_number_str, 
                     gene_id_str, transcript_id_str], 
                    [reading_frame, gene_name, exon_number, gene_id, 
                     transcript_id]):
                    # since new dic, add empty list, then append
                    ensdict[location][subkey].append(subvalue)
                count += 1
    f.close()
    # save pickle object to file.
    pkl = open(ensdictionary, 'wb')        #outputfile1 ='100_part_HS.pkl'            
    pickle.dump(ensdict, pkl)
    pkl.close()
    return ensdict, count

def openpkl(ensdictionary):
    ensdict = pickle.load( open(ensdictionary, 'rb') )
    return ensdict

def getgencode():
    gencode = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
        'ata':'I', 'atc':'I', 'att':'I', 'atg':'M',
        'aca':'T', 'acc':'T', 'acg':'T', 'act':'T',
        'aac':'N', 'aat':'N', 'aaa':'K', 'aag':'K',
        'agc':'S', 'agt':'S', 'aga':'R', 'agg':'R',
        'cta':'L', 'ctc':'L', 'ctg':'L', 'ctt':'L',
        'cca':'P', 'ccc':'P', 'ccg':'P', 'cct':'P',
        'cac':'H', 'cat':'H', 'caa':'Q', 'cag':'Q',
        'cga':'R', 'cgc':'R', 'cgg':'R', 'cgt':'R',
        'gta':'V', 'gtc':'V', 'gtg':'V', 'gtt':'V',
        'gca':'A', 'gcc':'A', 'gcg':'A', 'gct':'A',
        'gac':'D', 'gat':'D', 'gaa':'E', 'gag':'E',
        'gga':'G', 'ggc':'G', 'ggg':'G', 'ggt':'G',
        'tca':'S', 'tcc':'S', 'tcg':'S', 'tct':'S',
        'ttc':'F', 'ttt':'F', 'tta':'L', 'ttg':'L',
        'tac':'Y', 'tat':'Y', 'taa':'*', 'tag':'*',
        'tgc':'C', 'tgt':'C', 'tga':'*', 'tgg':'W'}
    return gencode
        
def translate(nucleotide_seq, reading_frame):
    gencode = getgencode()
    if reading_frame == '0':
        nucleotide_seq = nucleotide_seq
    elif reading_frame == '1':
        nucleotide_seq = nucleotide_seq[1:]
    elif reading_frame == '2':
        nucleotide_seq = nucleotide_seq[2:]
    elif reading_frame == '.':
        pass
    elif reading_frame == '':
        pass
    else:
        pass
    if reading_frame == '0' or reading_frame == '1' or reading_frame == '2':
        amino_acid_seq = ''.join([gencode.get(nucleotide_seq[3*i:3*i+3], 'X') \
                                    for i in range(len(nucleotide_seq)//3)])
        return amino_acid_seq
    else:
        print 'Reading frame must be 0, 1 or 2.'
        sys.exit()
        
def write_seq_to_file(myseq, mywritefile, maxlength=60):
    '''
    Takes a sequence, writes to writefile.
    Breaks down sequence to new line every 60 
    characters.
    Input:
        myseq: a long string, probably > 60
        mywritefile: writefile object from open()
    '''
    # init counters
    letter_count = 0
    seq_list = []
    for aa in myseq:
        seq_list.append(aa)
        letter_count = letter_count + 1
        if letter_count % maxlength == 0:
            # Max lenght reached, write file.
            mywritefile.write(''.join(seq_list + ['\n']))
            # reset counters
            seq_list = []
    # Done iterating, write remainders to file
    # only if seq_list is not empty.
    if len(seq_list) != 0:
        mywritefile.write(''.join(seq_list + ['\n']))
    return None

def extract_location_from_miso_event(miso_event, exon_number):
    '''
    Given MISO coordinate (e.g.:
    chr2:25650404:25650500:-@chr2:25642384:25642404:-@chr2:25611071:25611230:-
    Return either 1st, 2nd or 3rd (depending on exon_number) coordinate range.
    
    Exon number should be an integer.
    '''
    # split by @
    miso_coords_list = miso_event.split('@')
    
    # Get exon number, since index starts at 0, subtract exon_number by 1.
    return miso_coords_list[exon_number -1]

def main():
    ensdictionary = sys.argv[1]
    dna_fasta = sys.argv[2]
    try:
        exon_isoform_number = int(sys.argv[3])    # exon 1, 2 or 3...
    except ValueError:
        print '3rd argument must be integer, eg 1 2 or 3 (exon number)'
        sys.exit()
    summary_output_path = sys.argv[4]
    
    if len(sys.argv) < 5:
        print '4 arguments must be specified in command line.'
        sys.exit()
    
    print 'Loading dictionary from %s' %ensdictionary
    ensembl_dic = pickle.load(open(ensdictionary, "rb"))
    print 'Loaded dictionary from %s' %ensdictionary
    
    # init summary file, write header
    summary_file = open(summary_output_path, 'wb')
    summary_writer = csv.writer(summary_file, delimiter='\t')
    # 1. Get subkeys for ensembl dictionary
    reading_frame_str, \
    gene_name_str, \
    exon_number_str, \
    gene_id_str, \
    transcript_id_str, \
    _ = get_subkeys()
    # 2. Define additional column names in summaryfile
    miso_event_str = 'miso_event'
    nucleotide_seq_str = 'nucleotide_sequence'
    amino_acid_seq_str = 'amino_acid_sequence'
    # 3. write header
    outheader = [gene_name_str, miso_event_str, reading_frame_str, 
                nucleotide_seq_str, 
                 amino_acid_seq_str, gene_id_str, 
                 transcript_id_str,
                 exon_number_str]
    summary_writer.writerow(outheader)
    # end init summary file
    
    print 'Reading DNA fasta files and translating to protein...'
    writecount = 0
    with open(dna_fasta, 'rb') as dnafile:
        for line in dnafile:    
            # alternates header and sequence.
            # remove \n from end of each line
            line = line.strip()
            if line.startswith('>'):    # header of fasta...
                miso_event = line[1:]    # removes first '>' in header
                location = \
                    extract_location_from_miso_event(miso_event, 
                                                          exon_isoform_number)
                    
                # Get reading frames (as a list) and get all possible
                # translations.
                # nothing we can do if location not in dic...
                if location in ensembl_dic:
                    dna_seq = dnafile.next().strip()
                    # Get gene, transcript, exon annotations from dictionary
                    gene_names = \
                        list(set(ensembl_dic[location][gene_name_str]))
                    exon_numbers = \
                        ensembl_dic[location][exon_number_str]
                    gene_ids = \
                        list(set(ensembl_dic[location][gene_id_str]))
                    transcript_ids = \
                        list(set(ensembl_dic[location][transcript_id_str]))
                    
                    # Get non-redundant reading frames for iterating
                    reading_frames = \
                        list(set(ensembl_dic[location]['reading_frame']))
                    # print 'Gene name: %s' %ensembl_dic[location]['gene_name']
                    # print dna_seq
                    for reading_frame in reading_frames:
                        # print 'Reading frames: %s' %reading_frame
                        amino_acid_seq = translate(dna_seq, reading_frame)
                        # write to summary output path (order matches header)
                        # some are lists, so convert them to CSV
                        summary_writer.writerow([','.join(gene_names), 
                                                 miso_event, 
                                                 reading_frame,
                                                 dna_seq,
                                                 amino_acid_seq, 
                                                 ','.join(gene_ids), 
                                                 ','.join(transcript_ids),
                                                 ','.join(exon_numbers)])
                        writecount += 1
    print '%s rows written to: %s' %(writecount, summary_output_path)
                    
                    
if __name__ == '__main__':
    main()
