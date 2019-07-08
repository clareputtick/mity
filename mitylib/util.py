import logging
import subprocess

def tabix(f):
    """
    Generate a tabix index for a bgzipped file
    
    :param f: path to a bgzip compressed file
    :type f: str
    
    :returns: Nothing
    :rtype: None
    """
    tabix_call = "tabix " + f
    logging.debug(tabix_call)
    subprocess.run(tabix_call, shell=True)

