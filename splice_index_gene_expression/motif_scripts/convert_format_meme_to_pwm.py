'''
Created on 2013-09-12

@author: jyeung

Takes text file of meme format, converts to 
PWM format, so it is insertable into CISBP-RNA database.
'''


import sys
from motif_utils import motif_obj, convert_meme_to_pwm, write_pwm_obj_to_file, convert_chipmunk_to_pwm, convert_weeder_to_pwm


def main(motif_file, pwm_file):
    motif = motif_obj(motif_file)
    tool_used = motif.guess_motif_tool()
    if tool_used == 'Meme':
        pwm_obj = convert_meme_to_pwm(motif.file)
    elif tool_used == 'ChiPMunk':
        pwm_obj = convert_chipmunk_to_pwm(motif.reader)
    elif tool_used == 'Weeder':
        pwm_obj = convert_weeder_to_pwm(motif.reader)
    else:
        print('%s not recognized as Meme, ChiPMunk or Weeder.')
        sys.exit()
    write_pwm_obj_to_file(pwm_obj, pwm_file)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('motif file (input) and output file must be '\
              'specified in command line.')
        sys.exit()
    motif_file = sys.argv[1]
    pwm_file = sys.argv[2]
    main(motif_file, pwm_file)