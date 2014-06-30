'''
Created on 2014-01-24

@author: jyeung

Create motif file for homer.
Converts meme db to homer db.
Reads from output of build_custom_meme_db or any meme db
'''

import sys
from optparse import OptionParser
from motif_scripts.utilities.db_objs import get_seq_from_pwm,\
    write_motif_to_homer
from motif_scripts.utilities import db_objs


def main():
    usage = 'usage: %prog meme_db homer_db\n'\
        'Two args must be specified in commandline: \n'\
        '1) MEME database (input)\n'\
        '2) homer database (output)\n'\
        '-h to display this help message.\n'
    parser = OptionParser(usage=usage)
    parser.add_option('-s', '--log_score', dest='log_score',
                      help='Log Score for each motif, '\
                      'default is 7.5',
                      default='7.5')
    (options, args) = parser.parse_args()
    if len(args) != 2:
        print 'Incorrect number of parameters specified.'
        print usage
        sys.exit()
    input_db = args[0]
    output_db = args[1]
    log_score = options.log_score
    # check log-score is a number
    try:
        float(log_score)
    except ValueError:
        print 'Log score must be a float.'
    
    # init write file
    writefile = open(output_db, 'wb')
    
    with open(input_db, 'rb') as readfile:
        motif_obj = db_objs.MemeDB(readfile)
        n_motifs = 0
        while not motif_obj.end_of_file:    # itereate until end of file
            motif_id, pwm = motif_obj.get_next_motif()
            if motif_id is not None and pwm is not None:
                seq = get_seq_from_pwm(pwm)
                # write to file
                write_motif_to_homer(motif_id, seq, pwm, log_score, writefile)
                n_motifs += 1
    
    print '%s motifs created in: %s' %(n_motifs, output_db)
        
    # close
    writefile.close()

if __name__ == '__main__':
    main()