'''
Created on 2013-10-02

@author: jyeung
Requires selenium 2.0 and BeautifulSoup 4.3.1!
'''


# For BeautifulSoup test
import sys
import csv
import os
import urllib
import urllib2
import string
from BeautifulSoup import BeautifulSoup
from selenium import webdriver

def selenium_test():
    # Firefox Webdriver
    driver = webdriver.Firefox()    # no need to set things up like Chrome.
    
    # Go to codepad.orgse
    driver.get('http://codepad.org')
     
    # Select the Python language option
    python_link = driver.find_elements_by_xpath("//input[@name='lang' and @value='Python']")[0]
    python_link.click()
     
    # Enter some text!
    text_area = driver.find_element_by_id('textarea')
    text_area.send_keys("print 'Hello,' + ' World!'")
     
    # Submit the form!
    submit_button = driver.find_element_by_name('submit')
    submit_button.click()
     
    # Make this an actual test. Isn't Python beautiful?
    assert "Hello, World!" in driver.get_page_source()
     
    # Close the browser!
    driver.quit()

def beautiful_soup_test():
    user_agent = 'Mozilla/5 (Solaris 10) Gecko'
    headers = { 'User-Agent' : user_agent }
    values = {'s' : 'China' }
    data = urllib.urlencode(values)
    request = urllib2.Request("http://www.dict.cc/", data, headers)
    response = urllib2.urlopen(request)
    the_page = response.read()
    pool = BeautifulSoup(the_page)
    results = pool.findAll('td', attrs={'class' : 'td7nl'})
    source = ''
    translations = []
    
    for result in results:
        word = ''
        for tmp in result.findAll(text=True):
            word = word + " " + unicode(tmp).encode("utf-8")
        if source == '':
            source = word
        else:
            translations.append((source, word))
    
    for translation in translations:
        print "%s => %s" % (translation[0], translation[1])
        
def select_dropdown_menu(driver, menu_name, select_option):
    '''
    From menu name, select an option.
    '''
    # Create string for searching element.
    jstr = "//select[@name='%s']/option[text()='%s']" % (menu_name, 
                                                          select_option)
    driver.find_element_by_xpath(jstr).click()
    return jstr

def find_textbox(driver, textbox_id):
    input_element = driver.find_element_by_id(textbox_id)
    return input_element

def convert_tab_to_four_spaces(long_str):
    '''
    Convert str containing \t to four spaces.
    '''
    splittabs = long_str.split('\t')
    fourspaces = '    '.join(splittabs)
    return fourspaces

def read_textfile_as_string(textfile):
    '''
    Reads textfile and returns an object that is 
    easily insertable to selenium driver.
    '''
    with open(textfile, 'rU') as myfile:
        mylines = myfile.readlines()
        mylines_join = ''.join(mylines)
        mylines_fourspaces = convert_tab_to_four_spaces(mylines_join)
    return mylines_fourspaces

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

def retrieve_rbp_results(pagesource, tag, rbp_list):
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
    return list(set(rbp_results))

def write_list_to_file(list, output_filename):
    '''
    From a list, write to file, each list is its own row.
    '''
    with open(output_filename, 'wb') as myfile:
        jwriter = csv.writer(myfile, delimiter='\t')
        count = 0
        for i in list:
            jwriter.writerow([i])
            count += 1
    return count

def get_files_with_extension(folder, extension):
    '''
    Look in folder, get all folders with extension, (e.g. .pwm files)
    '''
    file_list = [f for f in os.listdir(folder) if f.endswith(extension)]
    return file_list

def query_cisrbp_get_rbp(driver, website, menuname, select_option, 
                         textbox_id, input_lines, tag, rbp_list):
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
    my_rbps = retrieve_rbp_results(driver.page_source, tag, rbp_list)
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
    
def main():
    # beautiful_soup_test()
    # selenium_test()
    pwm_folder = sys.argv[1]
    rbp_annotations = sys.argv[2]
    output_path = sys.argv[3]
    find_rbps_matching_motifs(pwm_folder, rbp_annotations, output_path)

if __name__ == '__main__':
    main()