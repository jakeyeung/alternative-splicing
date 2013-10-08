'''
Created on 2013-10-04

@author: jyeung

From fasta files, submit jobs to MEME.
'''

import os
import sys
from optparse import OptionParser
from selenium import webdriver
import webtool_utilities


def set_output_email(driver, email):
    '''
    meme requires an email so they can send you a link to the results.
    '''
    # init textbox names
    emailbox = 'email'
    confirm_emailbox = 'email_confirm'
    
    for txtbox in [emailbox, confirm_emailbox]:
        webtool_utilities.submit_text_to_textbox(driver, email, txtbox)
    return None

def set_meme_parameters(driver, min_width, max_width, n_motifs):
    '''
    Select meme paramters. Currently can set shortest word (min_width),
    longest word (max_width), number of motifs (n_motifs), single_strand or not
    
    MEME defaults are to be cleared.
    '''
    # Init textbox names
    minw_box = 'minw'
    maxw_box = 'maxw'
    nmotifs_box = 'nmotifs'
    
    for txt, txtbox in zip([min_width, max_width, n_motifs], 
                           [minw_box, maxw_box, nmotifs_box]):
        webtool_utilities.submit_text_to_textbox(driver, txt, txtbox, 
                                                 clear_defaults=True)
    return None

def set_meme_options(driver, dist, is_iss):
    '''
    Set distribution and whether to it is single stranded or not.
    
    Distribution options:
    one, zero_one, any.
    one: expect motif to show up once per sequence
    zero_one: expect motif to show up zero or one per sequence
    any: expect motif to show up any number of times per sequence
    Default to any if unsure.
    
    is_iss: Recommend True for RNA-motifs. False for DNA.
    '''
    # Init option names
    if dist == 'any':
        dist_val = 'anr'
    elif dist == 'one':
        dist_val = 'oops'
    elif dist == 'zero_one':
        dist_val = 'zoops'
    else:
        print('Distribution must be "any", "one", or '\
              '"zero_one". %s found.' %dist)
        sys.exit()
    is_iss_checkbox_name = 'posonly'
    
    # Click on radio button for distribution
    dist_radio_butt = \
        driver.find_element_by_xpath('//input[@value="%s"]' %dist_val)
    dist_radio_butt.click()
    
    # Checkbox for is_single_strand
    is_iss_checkbox = driver.find_element_by_name(is_iss_checkbox_name)
    is_iss_checkbox.click()
    return None

def set_meme_description(driver, descript):
    # Set box names
    descrip_name = 'description'
    descrip_box = driver.find_element_by_name(descrip_name)
    descrip_box.send_keys(descript)
    return None

def main(fasta_file_list, min_width, max_width, n_motifs, 
         email, dist, is_iss):
    # Set constants
    website = 'http://meme.nbcr.net/meme/cgi-bin/meme.cgi'
    textbox_id_name = 'data'
    submit_button_name = 'target_action'
    
    # Initialize web driver
    driver = webdriver.Firefox()
    
    for fasta_file in fasta_file_list:
        # Get fasta files for input to meme
        fasta_lines = webtool_utilities.read_textfile_as_string(fasta_file)
        print 'Number of characters: %s' %len(fasta_lines)
        if len(fasta_lines) > 70000:
            print('Fasta file: %s possibly too large.\n'\
                  'Warning: number of characters (%s) may be too large for meme to '\
                  'handle.\nDouble check that the fasta file '\
                  'was actually submitted to meme.' %(fasta_file, 
                                                      len(fasta_lines)))
        # Go to meme website
        driver.get(website)
        
        # Set email address to be submitted
        set_output_email(driver, email)
        
        # Set description to be fasta filename.
        descript = os.path.basename(fasta_file)
        set_meme_description(driver, descript)
        
        # Select meme settings
        set_meme_parameters(driver, min_width, max_width, n_motifs)
        
        # Select meme options
        set_meme_options(driver, dist, is_iss)
        
        # Submit fasta lines to website.
        webtool_utilities.submit_text_to_textbox(driver, fasta_lines, 
                                                 textbox_id_name)
        
        # Submit
        driver.find_element_by_name(submit_button_name).click()
    
    driver.close()
    
if __name__ == '__main__':
    
    usage = 'usage: %prog [options] fasta_file email'
    parser = OptionParser(usage=usage)
    
    parser.add_option('-s', '--min_width', dest='min_width', default=4,
                      help='Length of shortest motif to be discovered '\
                      '(Default 4).')
    parser.add_option('-l', '--max_width', dest='max_width', default=21,
                      help='Length of longest motif to be discovered '\
                      '(Default 22).')
    parser.add_option('-N', '--n_motifs', dest='n_motifs', default=10,
                      help='Number of motifs to be discovered, default 10')
    parser.add_option('--distribution', dest='dist', default='zero_one',
                      help='Occurence of a single motif is distributed how?'\
                      'Possible values: one, zero_one, any. Default: zero_one')
    parser.add_option('--is_single_strand', dest='is_ss', default=True,
                      help='Search given strand only? Default True.')
    parser.add_option('--batch_mode', dest='batch_mode', default=False,
                      help='If True: submits all fasta files in fasta_file directory to MEME.\n'\
                      'If False (Default): only submits fasta file in arg1 to MEME.\n'\
                      'Note: Keep arg1 as fastas_file, rather than directory even in batch_mode')
    (options, args) = parser.parse_args()
    
    fasta_file = args[0]
    email = args[1]
    
    if options.batch_mode==True:
        fasta_dir = os.path.dirname(fasta_file)
        fasta_file_list = \
            [os.path.join(fasta_dir, f) for f in os.listdir(fasta_dir) \
                     if f.endswith('.fasta')]
    else:
        fasta_file_list = [fasta_file]
    
    main(fasta_file_list, options.min_width, options.max_width, 
         options.n_motifs, email, options.dist, options.is_ss)
    



