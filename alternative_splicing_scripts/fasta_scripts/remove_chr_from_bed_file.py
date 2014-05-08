'''
Created on 2013-11-26

@author: jyeung
Bed files get column "chr1" to be "1".
'''

import sys
import csv
from optparse import OptionParser

def main():
    # Get parse options.
    usage = 'usage: %prog input_bed_file output_bed_file.\n'\
        'Following must be specified in command line:\n'\
              'input_bed_file: Bedfiles with chr-names like "chr1"'\
                '\n'\
              'output_bed_file: bedfile will have chr-names like "1"\n'
    parser = OptionParser(usage=usage)
    (_, args) = parser.parse_args()
    if len(args) < 2:
        print(usage)
        sys.exit()
    # Define args to variables
    input_bedfile = args[0]
    output_bedfile = args[1]
    
    # Init read and write file.
    writefile = open(output_bedfile, 'wb')
    mywriter = csv.writer(writefile, delimiter='\t')
    writecount = 0
    with open(input_bedfile, 'rb') as readfile:
        # write header
        writefile.write(readfile.readline())
        # loop through readfile, replace chr with '' in first column.
        myreader = csv.reader(readfile, delimiter='\t')
        for row in myreader:
            chromo = row[0]
            # Remove "chr" from chromosome name.
            row[0] = chromo.replace('chr', '')
            # Write new row to writefile
            mywriter.writerow(row)
            writecount += 1
    writefile.close()
    print '%s rows written to %s.' %(writecount, output_bedfile)

if __name__ == '__main__':
    main()