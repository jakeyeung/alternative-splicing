'''
Created on 2013-10-04

@author: jyeung

From fasta files, submit jobs to MEME.
'''

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
        submit_text_to_textbox(driver, email, txtbox)
    return None

def submit_text_to_textbox(driver, text, textbox, clear_defaults=False):
    '''
    1) Find textbox
    2) Insert text
    Clears defaults if set to true.
    '''
    input_element = driver.find_element_by_name(textbox)
    if clear_defaults:
        input_element.clear()
    input_element.send_keys(text)
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
        submit_text_to_textbox(driver, txt, txtbox, clear_defaults=True)
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
    dist_val = 'anr'
    is_iss_checkbox_name = 'posonly'
    
    # Click on radio button for distribution
    dist_radio_butt = \
        driver.find_element_by_xpath('//input[@value="%s"]' %dist_val)
    dist_radio_butt.click()
    
    # Checkbox for is_single_strand
    is_iss_checkbox = driver.find_element_by_name(is_iss_checkbox_name)
    is_iss_checkbox.click()
    
def submit_to_meme(driver):
    pass

def main(fasta_file, min_width, max_width, n_motifs, email, dist, is_iss):
    # Set constants
    website = 'http://meme.nbcr.net/meme/cgi-bin/meme.cgi'
    textbox_id_name = 'data'
    submit_button_name = 'target_action'
    
    # Initialize web driver
    driver = webdriver.Firefox()
    
    # Get fasta files for input to meme
    fasta_lines = webtool_utilities.read_textfile_as_string(fasta_file) 
    
    # Go to meme website
    driver.get(website)
    
    # Set email address to be submitted
    set_output_email(driver, email)
    
    # Select meme settings
    set_meme_parameters(driver, min_width, max_width, n_motifs)
    
    # Select meme options
    set_meme_options(driver, dist, is_iss)
    
    # Submit fasta lines to website.
    submit_text_to_textbox(driver, fasta_lines, textbox_id_name)
    
    # Submit
    driver.find_element_by_name(submit_button_name).click()
    
if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option('-f', '--file', dest='fasta_file', 
                      help='FASTA file from which to discover motifs')
    parser.add_option('-s', '--min_width', dest='min_width', default=4,
                      help='Length of shortest motif to be discovered '\
                      '(Default 4).')
    parser.add_option('-l', '--max_width', dest='max_width', default=22,
                      help='Length of longest motif to be discovered '\
                      '(Default 22).')
    parser.add_option('-N', '--n_motifs', dest='n_motifs', default=10,
                      help='Number of motifs to be discovered, default 10')
    parser.add_option('-e', '--email', dest='email',
                      help='Email to which MEME will send results.')
    parser.add_option('-d', '--distribution', dest='dist', default='any',
                      help='Occurence of a single motif is distributed how?'\
                      'Possible values: one, zero_one, any. Default: any')
    parser.add_option('--is_single_strand', dest='is_ss', default=True,
                      help='Search given strand only? Default True.')
    (options, args) = parser.parse_args()
    main(options.fasta_file, options.min_width, options.max_width, 
         options.n_motifs, options.email, options.dist, options.is_ss)
    



