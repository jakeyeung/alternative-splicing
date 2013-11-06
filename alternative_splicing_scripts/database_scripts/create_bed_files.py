'''
Created on 2013-06-20

@author: jyeung
From fasta annotation file, create bed file.
'''


import sys
import csv


# Set colname constants
chr_colname = 'Chromosome'
strand_colname = 'Strand'
# start_colname = 'Unit1_start_chr'
# end_colname = 'Unit1_end_chr'


class read_write_files(object):
    '''
    Read annotation file and write output files.
    '''
    
    def __init__(self, read_fullpath, write_fullpath):
        '''
        Constructor
        '''
        self.readfullpath = read_fullpath
        self.writefullpath = write_fullpath
        self.readrowcount = 0
        self.writerowcount = 0
        
    def __enter__(self):
        '''
        Open read and write paths as read and write objects.
        '''
        self.readfile = open(self.readfullpath, 'rb')
        self.reader = csv.reader(self.readfile, delimiter='\t')
        
        self.writefile = open(self.writefullpath, 'wb')
        self.writer = csv.writer(self.writefile, delimiter='\t')
        
    def __exit__(self, exittype, exitvalue, exittraceback):
        '''
        Close the file.
        '''
        self.readfile.close()
        self.writefile.close()
        
    def readnext(self):
        row = self.reader.next()
        self.readrowcount += 1
        return row
    
    def writenext(self, row):
        self.writer.writerow(row)
        self.writerowcount += 1
        return None

def create_bed_files_from_annotation(annot_fullpath, output_fullpath, 
                                     chr_colname, start_colname, end_colname, id_colname, strand_colname, 
                                     trackname, chromosome_list):
    '''
    From annotation file, create a bed file. 
    '''
    
    read_write_obj = read_write_files(annot_fullpath, output_fullpath)
    with read_write_obj:
        # Write header
        writeheader = 'track name=%s description=\"ALEXA_hs_65_37j (Human - build37/hg19)\" useScore=1\n' %trackname
        read_write_obj.writefile.write(writeheader)
        # Read header
        readheader = read_write_obj.readnext()
        while True:
            try:
                row = read_write_obj.readnext()
            except StopIteration:
                print('%s rows read, %s written, breaking loop.' %(read_write_obj.readrowcount, 
                                                                   read_write_obj.writerowcount))
                break
            chromosome = ''.join(['chr', row[readheader.index(chr_colname)]])
            if chromosome in chromosome_list:
                '''
                Filter out werid chromosome names.
                '''
                start = row[readheader.index(start_colname)]
                end = row[readheader.index(end_colname)]
                rowid = row[readheader.index(id_colname)]
                score = 999
                strand = int(row[readheader.index(strand_colname)])
                if strand < 0:
                    strand = '-'
                elif strand > 0:
                    strand = '+'
                else:
                    sys.exit('Strand neither greater or less than 0, %s' %strand)
                bed_row = [chromosome, start, end, rowid, score, strand]
                read_write_obj.writenext(bed_row)
    

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Fasta annotation file and output file full path must be specified.')
        sys.exit()
    
    # Define system arguments
    annot_fullpath = sys.argv[1]
    output_fullpath = sys.argv[2]
    trackname = sys.argv[3]
    id_colname = sys.argv[4]
    start_colname = sys.argv[5]
    end_colname = sys.argv[6]
    
    chromosome_list = [''.join(['chr', str(suffix)]) for suffix in range(1, 23)]
    chromosome_list.append('chrX')
    chromosome_list.append('chrY')
    
    # Read annotation file, write a bed file.
    create_bed_files_from_annotation(annot_fullpath, output_fullpath,
                                     chr_colname, start_colname, end_colname, id_colname, strand_colname,
                                     trackname, chromosome_list)