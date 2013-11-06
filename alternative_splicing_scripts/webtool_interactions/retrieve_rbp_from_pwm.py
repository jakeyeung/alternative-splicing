'''
Created on 2013-10-02

@author: jyeung
Requires selenium 2.0 and BeautifulSoup 4.3.1!
'''


# For BeautifulSoup test
import sys
import csv
import os
from selenium import webdriver
from webtool_utilities import select_dropdown_menu, find_textbox, \
    read_textfile_as_string, write_list_to_file, get_files_with_extension
try:
    from bs4 import BeautifulSoup    # linux
except ImportError:
    from BeautifulSoup import BeautifulSoup    # windows

def extract_rbps_from_annotations(rbp_annotations):
    '''
    Read RBP annotations downloaded from 
    homo sapiens http://cisbp-rna.ccbr.utoronto.ca/bulk.php
    Default filename RBP_Information.txt
    
    Assumes it is all homo sapiens, so extracts everything
    from RBP_Name column.
    '''
    # define constants
    rbp_colname = 'RBP_Name'
    
    rbp_list = []
    
    with open(rbp_annotations, 'rb') as myfile:
        myreader = csv.reader(myfile, delimiter='\t')
        header = myreader.next()
        for row in myreader:
            rbp_list.append(row[header.index(rbp_colname)])
    return rbp_list

def remove_spaces_from_results(result):
    '''
    For CISRBP, RBP name usually precedes with a long 
    whitespace, remove those for easier matching.
    '''
    try:
        result = result.string.replace(' ', '')
    except AttributeError:
        pass
    return result

def retrieve_rbp_results(pagesource, rbp_list):
    '''
    Parse a page source, looking for a particular tag
    (e.g. for CISRBP, rbp names are in tag 'a').
    Matches each possible result with an rbp name,
    returns list of results that contain rbp name.
    
    Function at the moment creates triplicates for each
    RBP, due to the fact that for each RBP, there are
    three links to it.
    
    I just set the results, not elegant, but it's ok.
    '''
    
    # Def constants
    tag = 'a'    # <a> contains our RBP results.
    
    # Use BeautifulSoup to parse html file
    results_source = BeautifulSoup(pagesource)
    all_results = results_source.findAll(tag)
    
    # Iterate results, matching to rbp_list
    rbp_results = []
    for r in all_results:
        try:
            r_stripped = r.string.strip()
        except AttributeError:
            pass
        if r_stripped in rbp_list:
            rbp_results.append(r_stripped)
    rbp_results = list(set(rbp_results))
    return rbp_results

def query_cisrbp_get_rbp(driver, input_lines, rbp_list):
    '''
    Inputs:
        Driver: firefox preferably.
        Input lines: PWM matrix containing A C G U header.
        rbp_list, list containing known RBPs for matching.
    '''
    # Set constants
    website = 'http://cisbp-rna.ccbr.utoronto.ca/TFTools.php'
    menuname = 'scanMotifSpec'
    select_option = 'Homo_sapiens'
    textbox_id = 'scanPWM'
    
    # Use Selenium to insert data and retrieve pagesource.
    # go to CISBP site
    driver.get(website)
    # Select homo sapiens
    select_dropdown_menu(driver, menuname, select_option)
    # Insert pwm_path_lines to textbox
    input_element = find_textbox(driver, textbox_id)
    # Insert keys
    input_element.send_keys(input_lines)
    input_element.submit()
    # Use BeautifulSoup to parse html file
    my_rbps = retrieve_rbp_results(driver.page_source, rbp_list)
    return my_rbps
            
def find_rbps_matching_motifs(pwm_folder, rbp_annotations, output_path):
    # Use both selenium and beautifulsoup to navigate page
    # and also parse HTML data.
    
    # Set constants
    website = 'http://cisbp-rna.ccbr.utoronto.ca/TFTools.php'
    menuname = 'scanMotifSpec'
    select_option = 'Homo_sapiens'
    textbox_id = 'scanPWM'
    tag = 'a'    # <a> contains our RBP results.
    pwm_ext = '.pwm'
    
    # Set webdriver
    driver = webdriver.Firefox()
    
    # Find PWM files from pwm_folder
    pwm_file_list = get_files_with_extension(pwm_folder, pwm_ext)
    
    # Iterate for each pwm file:
    all_rbps = []
    for pwm_f in pwm_file_list:
        # Create path
        pwm_path = os.path.join(pwm_folder, pwm_f)
        
        # Read textfiles
        rbp_list = extract_rbps_from_annotations(rbp_annotations)
        input_lines = read_textfile_as_string(pwm_path)
        
        my_rbps = \
        query_cisrbp_get_rbp(driver, website, menuname, select_option, 
                                       textbox_id, input_lines, tag, rbp_list)
        all_rbps += my_rbps
    
    all_rbps = list(set(all_rbps))
    # Write matched RBPs to file.
    print('RBPs found: %s.' %len(all_rbps))
    if len(all_rbps) > 0:
        counts = write_list_to_file(all_rbps, output_path)
        print('%s rows written to: %s' %(counts, output_path))
    
    # Close
    driver.close()
    
def main(pwm_folder, rbp_annotations, output_path):
    find_rbps_matching_motifs(pwm_folder, rbp_annotations, output_path)

if __name__ == '__main__':
    pwm_folder = sys.argv[1]
    rbp_annotations = sys.argv[2]
    output_path = sys.argv[3]
    main(pwm_folder, rbp_annotations, output_path)
    