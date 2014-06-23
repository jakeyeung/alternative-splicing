'''
Created on 2014-02-21

@author: jyeung
'''

class AnchorOutput(object):
    '''
    Representation of anchor output.
    '''


    def __init__(self, anchor_fname):
        '''
        Initialize anchor output filename.
        
        Reads entire file, stores it as object.
        '''
        with open(anchor_fname, 'rb') as readfile:
            self.lines = readfile.readlines()
    
    def binding_regions(self):
        '''
        Check if there are any predicted binding regions
        by searching for '#   None.\n'
        
        if not None:
            Parse lines between 
            '# Predicted binding regions\n'
            and
            '# Filtered regions\n'
        
        Grab start, end and length for each predicted binding region.
        
        Return information as a dic.
        
        Dictionary form:
        {binding_region_1: {start: 19, end: 52, length: 34}, ...}
        '''
        outdic = {}
        
        # Check there are binding regions
        try:
            # there should be "None" in output file under predicted
            # binding regions, so try to find it, if found, return None.
            self.lines.index('#   None.\n')
            return None
        except ValueError:
            # otherwise continue.
            pass
        
        # Get row start and end of predicted binding region information.
        # increase by 2 to predicted binding region and headers rows
        row_i = self.lines.index('# Predicted binding regions\n') + 2
        
        while True:
            '''
            Each row should be of form:
            #    1                        19          52          34
            Defined by:
            #    No.                   Start         End      Length
            contains mix of spaces and tabs. 
            Strategy: remove spaces, split by tab, return everything except
            first element (the #), also check they are integers.
            
            Break out of loop if row only contains #\n
            '''
            row = self.lines[row_i]
            # check it does not contain only a #
            if row == '#\n':
                break
            
            # remove \n, spaces, split by \t
            row = row.strip()
            row = row.replace(' ', '')
            row = row.split('\t')
            # first element is "#", ignore it.
            try:
                region_number = int(row[1])
                start = int(row[2])
                end = int(row[3])
                length = int(row[4])
            except ValueError:
                print 'Expected %s, %s, %s, %s to be integers.' \
                    %(row[1], row[2], row[3], row[4])
            # write to dic
            key = '_'.join(['binding_region', str(region_number)])
            # intialize subdic
            outdic[key] = {}
            for subkey, subval in zip(['start', 'end', 'length'], 
                                      [start, end, length]):
                outdic[key][subkey] = subval
            # increase row count and continue.
            row_i += 1
        return outdic
    
        