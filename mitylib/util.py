import logging
import subprocess

def tabix(outfile):
    """
    Generate a tabix index for a bgzipped file
    
    :param outfile: path to a bgzip compressed file
    :type outfile: str
    
    :returns: Nothing
    :rtype: None
    """
    tabix_call = "tabix " + outfile
    logging.debug(tabix_call)
    subprocess.run(tabix_call, shell=True)

