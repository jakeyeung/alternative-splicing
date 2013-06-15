'''
Created on 2013-06-14

@author: jyeung
This piece of code does several things:
1) Create more meaning headers for fasta files, converting headers called 1 to 
an id containing ensembl_name and also seq_name concatenated.
2) Increment alexa start position by one, deincrement alexa end position by one.
    Rationale: the raw start/end positions do not exactly match the ones in UCSC 
        unless we increment start position by one, deincrement end position by one.
3) Write these changes into a new annotation and fasta file. 
'''

import sys
import csv


class read_write_fasta_annot(object):
    '''
    For reading and writing fasta and annotation information.
    '''
    
    def __init__(self, fasta_read_fullpath, annot_read_fullpath, 
                 fasta_write_fullpath, annot_write_fullpath):
        '''
        Constructor
        '''
        self.fastareadpath = fasta_read_fullpath
        self.annotreadpath = annot_read_fullpath
        self.fastawritepath = fasta_write_fullpath
        self.annotwritepath = annot_write_fullpath
        self.fastareadrowcount = 0
        self.annotreadrowcount = 0
        self.fastawriterowcount = 0
        self.annotwriterowcount = 0
        
    def __enter__(self):
        '''
        Open read and write paths as read and write objects.
        '''
        self.fastareadfile = open(self.fastareadpath, 'rb')
        self.annotreadfile = open(self.annotreadpath, 'rb')
        self.fastareader = csv.reader(self.fastareadfile, delimiter='\t')
        self.annotreader = csv.reader(self.annotreadfile, delimiter='\t')
        
        self.fastawritefile = open(self.fastawritepath, 'wb')
        self.annotwritefile = open(self.annotwritepath, 'wb')
        self.fastawriter = csv.writer(self.fastawritefile, delimiter='\t')
        self.annotwriter = csv.writer(self.annotwritefile, delimiter='\t')
        
    def __exit__(self, exittype, exitvalue, exittraceback):
        '''
        Close the file.
        '''
        self.fastareadfile.close()
        self.annotreadfile.close()
        self.fastawritefile.close()
        self.annotwritefile.close()
    
    def readfastanext(self):
        row = self.fastareader.next()
        self.fastareadrowcount += 1
        return row
    
    def readannotnext(self):
        row = self.annotreader.next()
        self.annotreadrowcount += 1
        return row
        
    def writefastanext(self, row):
        self.fastawriter.writerow(row)
        self.fastawriterowcount += 1
        return None
    
    def writeannotnext(self, row):
        self.annotwriter.writerow(row)
        self.annotwriterowcount += 1
        return None

def create_ids(fasta_fullpath, annot_fullpath,
               fasta_output_fullpath, annot_output_fullpath):
    '''
    Read fasta and annot files line by line:
    
    In annotation file, read values from EnsEMBL_Gene_ID and Seq_Name
    columns, concatenate the two strings together, and replace first column
    (Junction_ID or Boundary_ID) with concatenated string (example, it would
    replace Junction_ID '1' with 'ENSG00000000003_E2a_E2d'. 
    
    In fasta file, replace header in each entry which is >1, >2, >3 ... 
    with the concatenated string (example, it would replace
    >1 with 'ENSG00000000003_E2a_E2d'
    '''
    # Set column names (in annot file) as constants.
    ensembl_id_colname = 'EnsEMBL_Gene_ID'
    seq_name_colname = 'Seq_Name'
    start_pos_colname = 'Unit1_start_chr'
    end_pos_colname = 'Unit1_end_chr'
    
    # Create read write obj to iterate and write rows in fasta and annot.
    read_write_obj = read_write_fasta_annot(fasta_fullpath, annot_fullpath, 
                                            fasta_output_fullpath, 
                                            annot_output_fullpath)
    with read_write_obj:
        # Get colnames of annotation file, write same colnames in annot output.
        annot_colnames = read_write_obj.readannotnext()
        read_write_obj.writeannotnext(annot_colnames)
        while True:
            # Itereate until no more rows to iterate in annotation file.
            try:
                annot_row = read_write_obj.readannotnext()
            except StopIteration:
                print('No more rows to iterate in annotation file.')
                print('%s and %s rows read in annot and fasta, '\
                      '%s and %s rows written in annot and fasta output.' \
                      %(read_write_obj.annotreadrowcount, 
                        read_write_obj.fastareadrowcount, 
                        read_write_obj.annotwriterowcount, 
                        read_write_obj.fastawriterowcount))
                break
            
            # Read ensembl_id, seq_name
            ensembl_id = annot_row[annot_colnames.index(ensembl_id_colname)]
            seq_name = annot_row[annot_colnames.index(seq_name_colname)]
            #
            # Create new id by concatenating the two strings
            #
            # Change boundary_ID or junction_ID (first column in annotrow) to
            # new_id
            new_id = '_'.join((ensembl_id, seq_name))
            annot_row[0] = new_id
            
            # increment start position by one, deincrement end pos by one.
            # then replace old start/end pos with new start/end pos.
            start_pos = int(annot_row[annot_colnames.index(start_pos_colname)])
            end_pos = int(annot_row[annot_colnames.index(end_pos_colname)])
            start_pos += 1
            end_pos -= 1
            annot_row[annot_colnames.index(start_pos_colname)] = start_pos
            annot_row[annot_colnames.index(end_pos_colname)] = end_pos
            
            # Write new annot_row to annot output
            read_write_obj.writeannotnext(annot_row)
            
            # Change fasta header to new_id, then write new header 
            # and its corresponding sequence to fasta output.
            try:
                fasta_header = read_write_obj.readfastanext()
                fasta_sequence = read_write_obj.readfastanext()
            except StopIteration:
                sys.exit('Ran out of rows to iterate in fasta file.')
            # Check first character in fasta_header is a '>'
            # If so, then proceed by replacing fasta_header with new_id.
            if fasta_header[0][0] == '>':
                fasta_header = ''.join(('>', new_id))
            else:
                sys.exit('fasta header does not begin with >, %s found.' %fasta_header[0][0])
            # Write header
            read_write_obj.writefastanext([fasta_header])
            # Write sequence (line after header)
            read_write_obj.writefastanext(fasta_sequence)
            

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('exonBoundaries fasta, annotations,  '\
              'and output filename must be '\
              'specified in commandline.')
        sys.exit()
    
    # Define system arguments
    fasta_fullpath = sys.argv[1]
    annot_fullpath = sys.argv[2]
    fasta_output_fullpath = sys.argv[3]
    annot_output_fullpath = sys.argv[4]
    
    # Reannotate and add better IDs for annotation file and fasta file. 
    create_ids(fasta_fullpath, annot_fullpath,
               fasta_output_fullpath, annot_output_fullpath)
    
    
    
    