'''
Created on 2013-06-21

@author: jyeung

Find out the UNIQUE modules between two optdis outputs. 
'''


import sys
import csv


def get_modules_from_file(file_fullpath, number_of_modules):
    '''
    The function assumes the file_fullpath has no headers,
    first column contains gene names, second column contains
    a subnetwork 'number' which ranges from 1 to number_of_modules.
    
    Returns a list of lists containing modules
    '''
    done = False    # For looping.
    all_modules = []    # for appending
    current_gene_module = []    # Initialize for while loop
    current_mod_numb = 1    # Initialize for while loop.
    module_number_list = range(1, number_of_modules + 1)    # For looping.
    with open(file_fullpath, 'rb') as modulefile:
        modulereader = csv.reader(modulefile, delimiter='\t')
        for mod_numb in module_number_list:
            while current_mod_numb == mod_numb:
                try:
                    row = modulereader.next()
                except StopIteration:
                    print('%s modules iterated, breaking.' %mod_numb)
                    all_modules.append(current_gene_module)
                    done = True
                    break
                current_mod_numb = int(row[1])    # module number
                current_gene_module.append(row[0])    # gene name
            '''
            # Broken out of while loop, means we have iterated
            # to the next module number, so add current_gene_module
            # to modules list, but the last element in current_gene_module
            # actually belongs to the next module number, so 
            # remember to put that last element as first member of 
            # the next module. 
            '''
            if done == False:
                all_modules.append(current_gene_module[:-1])
                # Reinitialize current_gene_module, keeping last element. 
                current_gene_module = [current_gene_module[-1]]
            elif done == True:
                pass
            else:
                sys.exit('Done must be False or True.')
    return all_modules

def find_unique_and_common_modules(reference_list, target_list):
    '''
    Given a reference list, check a list of interest for unique modules, then
    return those unique and common modules. 
    '''
    common_modules_list = []
    unique_modules_list = []
    for module in target_list:
        if module in reference_list:
            common_modules_list.append(module)
        else:
            unique_modules_list.append(module)
    return unique_modules_list, common_modules_list

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('Two optdis output fullpaths and number of modules must be provided in commandline.')
        sys.exit()
    optdis_output1 = sys.argv[1]
    optdis_output2 = sys.argv[2]
    try:
        number_of_modules = int(sys.argv[3])
    except ValueError:
        sys.exit('Number of modules must be integer, not %s' %number_of_modules)
    
    # Get modules from file.
    modules1 = get_modules_from_file(optdis_output1, number_of_modules)
    modules2 = get_modules_from_file(optdis_output2, number_of_modules)
    
    # Find unique modules for modules 1 and 2.
    unique_modules1, common_modules1 = find_unique_and_common_modules(modules2, modules1)
    unique_modules2, common_modules2 = find_unique_and_common_modules(modules1, modules2)
    
    print len(unique_modules1), unique_modules1
    print len(common_modules1), common_modules1
    print len(unique_modules2), unique_modules2
    print len(common_modules2), common_modules2
    