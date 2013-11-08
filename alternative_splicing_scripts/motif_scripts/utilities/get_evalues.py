'''
Created on 2013-11-06

@author: jyeung
'''

import re

def retrieve_evalues_from_file(meme_file):
    matchcount = 0
    evalues = []
    with open(meme_file, 'rb') as readfile:
        for line in readfile:
            m = re.search('(?<=E-value = )\S+', line)
            if m:
                matchcount += 1
                evalues.append(m.group(0))
    if matchcount != 10:
        print 'Warning: expected 10 (%s found) Evalues to be retrieved from %s'\
             %(matchcount, meme_file)
    return evalues
