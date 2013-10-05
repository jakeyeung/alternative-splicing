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
    # init constants
    textbox_name = 'email'
    # Init textbox
    input_element = driver.find_element_by_name(textbox_name)
    # Insert keys
    input_element.send_keys(email)
    input_element.submit()

def submit_text_to_textbox(driver, text, textbox):
    '''
    Initialized a web driver, use that to go to website, submit
    input into textbox.
    '''
    pass

def set_meme_parameters(driver, min_width, max_width, 
                        n_motifs, single_strand=True):
    '''
    Select meme paramters. Currently can set shortest word (min_width),
    longest word (max_width), number of motifs (n_motifs), single_strand or not
    '''
    pass

def submit_to_meme(driver):
    pass

def main(fasta_file, min_width, max_width, n_motifs, email):
    # Set constants
    website = 'http://meme.nbcr.net/meme/cgi-bin/meme.cgi'
    select_option = 'Homo_sapiens'
    textbox_id = 'data'
    
    # Initialize web driver
    driver = webdriver.Firefox()
    
    # Get fasta files for input to meme
    fasta_lines = webtool_utilities.read_textfile_as_string(fasta_file)
    
    # Go to meme website
    driver.get(website)
    
    # Set email address to be submitted
    set_output_email(driver, email)
    
    # Submit fasta lines to website.
    submit_text_to_textbox(driver, fasta_lines, textbox_id)
    
    # Select meme settings
    set_meme_parameters(driver, min_width, max_width, n_motifs)
    
    # Submit
    submit_to_meme(driver)
    
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-f', '--file', dest='fasta_file', 
                      help='FASTA file from which to discover motifs')
    parser.add_option('-minw', '--min_width', dest='min_width',
                      help='Length of shortest motif to be discovered '\
                      '(recommend 4 for CISBP-RNA db).')
    parser.add_option('-maxw', '--max_width', dest='max_width',
                      help='Length of longest motif to be discovered '\
                      '(recommend 22 for CISBp-RNA db).')
    parser.add_option('-N', '--n_motifs', dest='n_motifs',
                      help='Number of motifs to be discovered. Recommend 10')
    parser.add_option('-e', '--email', dest='email',
                      help='Email to which MEME will send results.')
    (options, args) = parser.parse_args()
    if len(args) != 5:
        parser.error('Must use five arguments -h or --help for tips.')
        sys.exit()
    main(options.fasta_file, options.min_width, options.max_width, 
         options.n_motifs, options.email)
    



