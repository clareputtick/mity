import logging
import subprocess
import os

def tabix(f):
    """
    Generate a tabix index for a bgzipped file
    
    :param f: path to a bgzip compressed file
    :type f: str
    
    :returns: Nothing
    :rtype: None
    """
    tabix_call = "tabix -f " + f
    logging.debug(tabix_call)
    subprocess.run(tabix_call, shell=True)

def check_missing_file(file_list, die=True):
    missing_files = []
    for item in file_list:
        if not os.path.isfile(item):
            missing_files.append(item)
    if die and len(missing_files) > 0:
        raise ValueError(f"Missing these files: {missing_files}")
    return missing_files

