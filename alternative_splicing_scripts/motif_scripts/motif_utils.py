'''
Created on 2013-09-12

@author: jyeung
'''

import csv 
import os

class stamp_obj(object):
    '''
    Object to load stamp file
    '''
    
    def __init__(self, stamp_file):
        '''
        Initialize reader
        '''
        self.file = open(stamp_file, 'rbU')
        
def remove_slash_n(name):
    '''
    Converts DE Meme1\n to DE Meme
    '''
    return name.split('\n')[0]
        
def create_filename_prefix(fname_prefix, name):
    '''
    Create a name like fname_prefix_DE_Meme1.txt
    name -> expects something like DE Meme1\n
    '''
    # There is \n, split by \n grab the first element.
    name = remove_slash_n(name)
    name_join = '_'.join(name.split(' '))
    filename = '.'.join([fname_prefix, name_join, 'txt'])
    return filename

def is_this_first_line(line):
    '''
    Expects first line to contain DE.
    '''
    if 'DE' in line:
        return True
    else:
        return False
    
def is_this_last_line(line):
    '''
    Expects last line to contain XX
    '''
    if 'XX' in line:
        return True
    else:
        return False

def write_files(myfile, fname_prefix, mydir):
    '''
    For every motif, write its own separate file.
    Make its name to be
    '''
    first_line = True    # intialize.
    for line in myfile:
        last_line = is_this_last_line(line)
        '''
        print line
        print('First line:%s' %first_line)
        print ('Last line:%s' %last_line)
        raw_input()
        '''
        if first_line:
            # First line: create file and write first row.
            filename = create_filename_prefix(fname_prefix, line)
            filename_path = os.path.join(mydir, filename)
            writefile = open(filename_path, 'wb')
            writefile.write(line)
            first_line = False
        elif last_line:
            writefile.close()
            print('File written to: %s' %filename_path)
            first_line = True
        else:
            writefile.write(line)

class motif_obj(object):
    '''
    Object to play with a motif. Convert to pwm etc.
    '''

    def __init__(self, motif_file):
        '''
        Constructor
        '''
        self.file = open(motif_file, 'rbU')
        self.header = self.file.readline()
        self.reader = csv.reader(self.file, delimiter=' ')
        
    def guess_motif_tool(self):
        '''
        Read the header and guess if motif is made from tools:
            meme
            weeder
            chipmunk
        Expects header to be a one liner, no other columns.
        '''
        # def constants
        meme_str = 'Meme'
        weeder_str = 'Weeder'
        chipmunk_str = 'ChiPMunk'
        for j_str in [meme_str, weeder_str, chipmunk_str]:
            if j_str in self.header:    # So it is unlisted. 
                print('%s is from %s tool.' %(self.header, j_str))
                self.tool = j_str
                break
        else:
            print('%s: no tool found.' %self.header)
            self.tool = None
        return self.tool

def write_pwm_header():
    '''
    Write first line of pwm format, which is A C G U.
    '''
    header = ['', 'A', 'C', 'G', 'U']
    return header
    
def convert_meme_to_pwm(readfile, header=False):
    '''
    Iterate rows in readobj, assuming it is meme format
    then make it into pwm format, which contains
    A C G U as first line and no bases on the last column.
    
    Store in list of lists which can be easily written to file.
    
    The darned thing is DOUBLE SPACE delimited! So the csv reader
    doesnt read the space delimited very well...
    '''
    # Init output list of lists
    pwm_list = []
    # Write header
    pwm_list.append(write_pwm_header())
    # Iterate rows
    for line in readfile:
        # Remove last column of every meme row.
        linesplit = line.split('  ')
        pwm_list.append(linesplit[0:-1])
    return pwm_list

def convert_chipmunk_to_pwm(reader, header=False):
    '''
    Converts chipmunk format to pwm format.
    
    To convert: divide each element by sum of the row
    to get probability.
    '''
    # Init output list of lists
    pwm_list = []
    # Write header
    pwm_list.append(write_pwm_header())
    # iterate rows
    for row in reader:
        # Last row is an empty space, ignore it.
        # First element is a counter, ignroe it.
        # Get weighted probability.
        count = row[0]
        row_floats = [float(i) for i in row[1:-1]]
        row_prob = [i/sum(row_floats) for i in row_floats]
        pwm_list.append([count] + row_prob)
    return pwm_list
        
def convert_weeder_to_pwm(reader, header=False):
    '''
    Format is the same as for chipmunk...
    '''
    pwm_list = convert_chipmunk_to_pwm(reader, header=False)
    return pwm_list
        
def write_pwm_obj_to_file(pwm_obj, output_file):
    '''
    Get pwm_obj from convert_x_to_pwm, then write
    that list of list to file.
    '''
    writecount = 0
    with open(output_file, 'wb') as outfile:
        outwriter = csv.writer(outfile, delimiter='\t')
        for row in pwm_obj:
            outwriter.writerow(row)
            writecount += 1
    print('%s rows written to file: %s' %(writecount, output_file))
    outfile.close()
    return writecount

def create_pwm_filename(motif_file_name):
    '''
    Take motif file name and append .pwm before the end.
    myfile.txt -> myfile.pwm.txt
    '''
    pwm_str = 'pwm'
    # Split
    motif_file_split = motif_file_name.split('.')
    # Insert
    motif_file_split.insert(-1, pwm_str)
    # Rejoin
    return '.'.join(motif_file_split[0:-1])        