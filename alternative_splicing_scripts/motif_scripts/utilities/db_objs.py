'''
Created on 2014-01-24

@author: jyeung
'''

def get_nuc_ind_dic():
    '''
    We want:
    Convert index (0, 1, 2, 3)
    to A C G T.
    
    This function makes a dic to easily do that.
    '''
    nuc_ind_dic = {0: 'A',
                   1: 'C',
                   2: 'G',
                   3: 'T'}
    return nuc_ind_dic

def get_seq_from_pwm(pwm):
    '''
    Given PWM, get the consensus sequence
    (i.e. get letters corresponding to highest
    probability).
    If both letters are equal probability, then choose
    alphabetically.
    
    Assumes element in PWM is in alphabetical order (A C G T)
    '''
    # nuc_ind_dic correlates index to nucleotide.
    nuc_ind_dic = get_nuc_ind_dic()
    
    # convert pwm to list of lists by splitting by tab
    pwm = [i.split('\t') for i in pwm]
    # convert to floats, but preserve list of list structure
    pwm = [[float(i) for i in sublist] for sublist in pwm]
    
    # init list, append letters to it corresopnding to nucleotide seq
    nuc_seq = []
    
    # get nucleotide corresopnding to highest fraction in each position.
    for row in pwm:
        # get index of largest value.
        i = row.index(max(row))
        nuc_seq.append(nuc_ind_dic[i])
    
    # Join list to string, return
    return ''.join(nuc_seq)
    
def write_motif_to_homer(motif_id, seq, pwm, log_score, writefile):
    '''
    Inputs:
        motif_id: e.g. CELF3,M004_0.6,I
        seq: nucleotide sequence representing largest probability at each pos
        pwm: list, each element is a row ready to be written
            to outfile.
        log_score: user defined, should be a number but string.
        writefile: write file object
    Outputs:
        None
        
    First line should be:
    >seq\tmotif_id\tlog_score
    Then pwm to follow.
    '''
    # add ">" to begining of seq
    seq = ''.join(['>', seq])
    # add "\n" to end of log_score
    log_score = ''.join([log_score, '\n'])
    # write header
    header = '\t'.join([seq, motif_id, log_score])
    writefile.write(header)
    
    # write pwm
    for line in pwm:
        writefile.write(line)
    
    return None

def get_motif_id_from_line(line):
    '''
    Given:
    MOTIF CELF3,M004_0.6,I
    Retrieve CELF3,M004_0.6,I
    
    Space delimited.
    
    need to strip it in case it has \n
    '''
    return line.strip().split(' ')[1]

class MemeDB(object):
    '''
    MEME file with a number of motifs.
    '''


    def __init__(self, readfile):
        '''
        Constructor
        '''
        self.file = readfile
        self.end_of_file = False
    
    def get_next_motif(self):
        '''
        Reads next motif by searching for
        line beginning with:
        MOTIF 
        then getting the PWM until blank line is retrieved.
        Return list of lists
        [[line1], [line2], ....]
        '''
        # init output which is list of lists.
        output_list = []
        for line in self.file:
            if line.startswith('MOTIF'):    # beginning of motif
                motif_line = line
                motif_id = get_motif_id_from_line(motif_line)
                # skip next 2 lines to begin PWM
                for _ in range(2):
                    motif_line = self.file.next()
                while not motif_line == '\n':
                    output_list.append(motif_line)
                    try:
                        motif_line = self.file.next()
                    except StopIteration:
                        break
                return motif_id, output_list
        # if it reaches here, it probably means it hit end of file, print
        # saying we reached end of file, return None, None and set 
        # end of file to True
        self.end_of_file = True
        return None, None
        