'''
Created on 2013-10-04

@author: jyeung
'''

import csv
import os
import urllib
import urllib2
from selenium import webdriver
try:
    from bs4 import BeautifulSoup    # linux
except ImportError:
    from BeautifulSoup import BeautifulSoup    # windows

def str_to_bool(true_false_str):
    '''
    Try to convert str to boolean. Works if string looks like True or False.
    '''
    if true_false_str in ['True', 'TRUE',' true', 't', 'T']:
        return True
    elif true_false_str in ['False', 'FALSE', 'false', 'f', 'F']:
        return False
    else:
        print('Option must be True or False: %s found.' %true_false_str)
        return None

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

def read_textfile_as_string(textfile, convert_tab_to_space=True):
    '''
    Reads textfile and returns an object that is 
    easily insertable to selenium driver.
    '''
    with open(textfile, 'rU') as myfile:
        mylines = myfile.readlines()
        mylines_join = ''.join(mylines)
        if convert_tab_to_space:
            mylines_join = convert_tab_to_four_spaces(mylines_join)
    return mylines_join

def write_list_to_file(mylist, output_filename):
    '''
    From a mylist, write to file, each mylist is its own row.
    '''
    with open(output_filename, 'wb') as myfile:
        jwriter = csv.writer(myfile, delimiter='\t')
        count = 0
        for i in mylist:
            jwriter.writerow([i])
            count += 1
    return count

def get_files_with_extension(folder, extension):
    '''
    Look in folder, get all folders with extension, (e.g. .pwm files)
    '''
    file_list = [f for f in os.listdir(folder) if f.endswith(extension)]
    return file_list

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