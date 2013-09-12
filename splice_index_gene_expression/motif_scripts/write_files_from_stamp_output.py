'''
Created on 2013-09-12

@author: jyeung
From stamp motif output, write a single text file
for every single motif.
'''

import sys
import os
from motif_utils import stamp_obj, write_files


def main(stamp_file, motif_name):
    stamp = stamp_obj(stamp_file)
    mydir = os.path.dirname(stamp_file)
    write_files(stamp.file, motif_name, mydir)

if __name__ == '__main__':
    if len(sys.argv) < 1:
        print('Stamp output file must be specified in command line.')
        sys.exit()
    stamp_file = sys.argv[1]
    motif_name = sys.argv[2]    # e.g. mincount_4_intron2_in
    main(stamp_file, motif_name)