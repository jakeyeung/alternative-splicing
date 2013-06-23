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
                current_gene_module.append((row[0], current_mod_numb))    # gene name
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
    Assumes reference list is a list of tuples, with first element as gene name
    second element is the module number.
    '''
    common_modules_list = []
    unique_modules_list = []
    
    # Remove the module number from tuple list, getting only a list of lists
    # with gene names to compare target_list gene modules to ref gene modules.
    # Why? Because you may have different gene module number but same genes.
    ref_gene_list = []
    for ref_tuplelist in reference_list:
        ref_gene_list.append([i[0] for i in ref_tuplelist])
        
    for tuplelist in target_list:
        if [i[0] for i in tuplelist] in ref_gene_list:
            common_modules_list.append(tuplelist)
        else:
            unique_modules_list.append(tuplelist)
    return unique_modules_list, common_modules_list

def write_modules_to_file(modules, output_fullpath):
    '''
    # Write unique modules first, then common modules.
    # Use same output_fullpath but add suffix in end to indicate
    # unique or common modules.
    '''
    with open(output_fullpath, 'wb') as writefile:
        writer = csv.writer(writefile, delimiter='\t')
        for tuple_list in modules:
            for tup in tuple_list:
                writer.writerow(tup)
    print('File written to %s ' %output_fullpath)
    return None

def find_super_unique_modules(modules, reference_modules):
    '''
    Super unique: modules where all of its members are unique to a 
    reference list.
    '''
    super_unique_modules = []    # initialize
    reference_list = []    # Flattened list of only genes
    for tup_list in reference_modules:
        for tup in tup_list:
            reference_list.append(tup[0])
    reference_list = set(reference_list)
    for tup_list in modules:
        # Find intersection between gene in tup_list and ref_list
        gene_list = [i[0] for i in tup_list]
        gene_ref_intersect = set(gene_list).intersection(reference_list)
        if len(gene_ref_intersect) == 0:
            super_unique_modules.append(tup_list)
    return super_unique_modules

if __name__ == '__main__':
    if len(sys.argv) < 6:
        print('Two optdis output fullpaths (for reading), number of modules, '\
              'and two write output paths (for writing) must be provided '\
              'in commandline.')
        sys.exit()
    optdis_output1 = sys.argv[1]
    optdis_output2 = sys.argv[2]
    try:
        number_of_modules = int(sys.argv[3])
    except ValueError:
        sys.exit('Number of modules must be integer, not %s' %number_of_modules)
    write_output1_unique = sys.argv[4]
    write_output1_common = sys.argv[5]
    write_output1_superunique = sys.argv[6]
    write_output2_unique = sys.argv[7]
    write_output2_common = sys.argv[8]
    write_output2_superunique = sys.argv[9]
    
    # Get modules from file.
    modules1 = get_modules_from_file(optdis_output1, number_of_modules)
    modules2 = get_modules_from_file(optdis_output2, number_of_modules)
    
    # Find unique modules for modules 1 and 2.
    unique_modules1, common_modules1 = find_unique_and_common_modules(modules2, modules1)
    unique_modules2, common_modules2 = find_unique_and_common_modules(modules1, modules2)
    
    # Are there TRULY unique modules?
    super_uniques2 = find_super_unique_modules(unique_modules2, modules1)
    super_uniques1 = find_super_unique_modules(unique_modules1, modules2)
    
    # Write unique and common modules to file. 
    for mod_tup_list, output_fullpath in zip([unique_modules1, common_modules1, 
                                              super_uniques1, 
                                          unique_modules2, common_modules2, 
                                          super_uniques2],
                                            [write_output1_unique,
                                             write_output1_common,
                                             write_output1_superunique,
                                             write_output2_unique,
                                             write_output2_common,
                                             write_output2_superunique]):
        write_modules_to_file(mod_tup_list, output_fullpath)
    