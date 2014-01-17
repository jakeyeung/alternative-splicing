'''
Created on 2014-01-17

@author: jyeung
'''

import sys


def get_inclusion_exclusion(incl_file=None, excl_file=None):
    '''
    Read inclusion and exclusion fasta file containing miso events,
    return dictionary containing 
    
    {miso_event: 'inclusion'|'exclusion' ,...}
    '''
    # Define keys
    incl_str = 'inclusion'
    excl_str = 'exclusion'
    outdic = {}
    
    for jfile, incl_or_excl in zip([incl_file, excl_file], 
                                   [incl_str, excl_str]):
        with open(jfile, 'rb') as readfile:
            for line in readfile:
                if line.startswith('>'):    # contains miso event
                    # remove first character '>', strip \n
                    miso_event = line[1:].strip()
                    if miso_event not in outdic:
                        outdic[miso_event] = incl_or_excl
                    else:
                        print 'Duplicate in miso event: %s' %miso_event
                        print 'Exiting...'
                        sys.exit()
    return outdic