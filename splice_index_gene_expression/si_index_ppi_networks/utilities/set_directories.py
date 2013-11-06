'''
Created on 2013-05-22

@author: jyeung
'''


import os


class mydirs(object):
    '''
    A class representing input, root and output folders. 
    '''
    def __init__(self, inputdir, outputdir):
        self.utilities = os.path.dirname(__file__)
        self.main = os.path.dirname(os.path.dirname(__file__))
        self.project = os.path.dirname(os.path.dirname(os.path.dirname(self.main)))
        self.inputdir = os.path.join(self.project, inputdir)
        self.outputdir = os.path.join(self.project, outputdir)
        
    
def set_directories(input_folder_name, cohort_folder_name, 
                    output_folder_name):
    '''
    Return current directory of MAIN (not utilities)
    Get project root dir
    Get input dir
    Get output dir
    Get plot dir
    '''
    _cur_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    _proj_dir = os.path.dirname(os.path.dirname(os.path.dirname(_cur_dir)))
    _input_dir = os.path.join(_proj_dir, input_folder_name, cohort_folder_name)
    _output_dir = os.path.join(_proj_dir, output_folder_name)
    
    return _cur_dir, _proj_dir, _input_dir, _output_dir