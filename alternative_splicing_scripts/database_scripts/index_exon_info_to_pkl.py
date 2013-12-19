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
    chr_list = [range(1, 23)] + ['X', 'Y']
    
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
    pkl = open(ensdictionary, 'wb')        #outputfile1 ='100_part_HS.pkl'
    #    Grabs readframe, constructs chromosome location, constructs exon 'identifiers' (transcriptnumber, jtype(like "exon"), number).
    
    # Define index names
    chr_index = 0    # chromosome, we want 1, 2, 3.. X and Y's only
    feature_index = 2    # we want coding sequence only
    exon_start_index = 3
    exon_end_index = 4
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
    for row in jreader:
        chr_name = row[chr_index]
        if chr_name not in chr_list:
            continue
        else:
            feature = row[feature_index]
            if feature != 'CDS':
                continue
            else:
                # Build key with chr:start:end
                start = row[exon_start_index]
                end = row[exon_end_index]
                location = ':'.join([chr_name, start, end])
                
                # get reading frame
                reading_frame = row[reading_frame_index]
                
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
                '''
                exon_id = \
                    extract_info_from_exon(exon_info, search=exon_id_str)
                '''
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
    f.close()            
    pickle.dump(ensdict, pkl)
    pkl.close()

def main():
    ensdatabase = sys.argv[1]
    ensdictionary = sys.argv[2]
    
    if len(sys.argv) < 3:
        print '3 arguments must be specified in command line.'
        sys.exit()
    
    print 'Making dictionary from %s' %ensdatabase
    make_dictionary(ensdatabase, ensdictionary)
    print 'Dictionary saved to %s' %ensdictionary
    
if __name__ == '__main__':
    main()
