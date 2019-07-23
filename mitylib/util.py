import logging
import subprocess
import os
import sys
import tempfile
import vcf
import configparser

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

def tmp_mity_file_name():
    """
    Create a tmp mity vcf file
    
    #TODO There must be more pythonic ways of doing this.
    """
    f = tempfile.NamedTemporaryFile(mode="wt", prefix='mity', suffix=".vcf",
                                    delete=False)
    f.close()
    assert isinstance(f.name, str)
    return f.name


def create_prefix(file_name, prefix=None):
    """
    Most mity functions have an optional prefix. If a prefix is not specified,
    then use the  file name (minus the .vcf.gz or .bam suffix) as the prefix.
    :param file_name: The vcf, bam, bed filename
    :param prefix: The optional prefix. If None, then craete a prefix from 
    vcf_name, else return prefix
    :return: str prefix
    """
    if ".vcf.gz" in file_name:
        prefix = [os.path.basename(file_name).split(".vcf")[0], prefix][prefix is not None]
    elif ".bam" in file_name:
        prefix = [os.path.basename(file_name).split(".bam")[0], prefix][prefix is not None]
    else:
        raise ValueError("Unsupported file type")
    return prefix


def write_vcf(new_vcf, out_file, genome_file='annot/b37d5.genome'):
    """
    write a vcf object to vcf.gz file with tbi index.
    
    This differs from write_merged_vcf, as mity merge doesn't split each VCF 
    line on tabs, whereas mity normalise does
    
    :param new_vcf: new_vcf is a list of lists, created by normalise
    :param out_file: the resulting filename. this should end in vcf.gz
    :return: None. This function writes a vcf.gz and vcf.gz.tbi file.
    """
    f = tmp_mity_file_name()
    logging.debug(f"Writing uncompressed vcf to {f}")
    with open(f, mode='wt', encoding='utf-8') as myfile:
        for vcf_line in new_vcf:
            myfile.write('\t'.join([str(elem) for elem in vcf_line]) + '\n')
    logging.debug(f"Sorting, bgzipping {f} -> {out_file}")
    subprocess.run(f"gsort {f} {genome_file} | bgzip -cf > {out_file}", shell=True)
    logging.debug(f"Tabix indexing {out_file}")
    tabix(out_file)
    os.remove(f)

def write_merged_vcf(new_vcf, out_file, genome_file='annot/b37d5.genome'):
    """
    write a vcf object to vcf.gz file with tbi index.
    
    This differs from write_vcf, as mity merge doesn't split each VCF line on 
    tabs, whereas mity normalise does
    
    :param new_vcf: new_vcf is a list of strings, created by merge
    :param out_file: the resulting filename. this should end in vcf.gz
    :return: None. This function writes a vcf.gz and vcf.gz.tbi file.
    """
    f = tmp_mity_file_name()
    logging.debug(f"Writing uncompressed vcf to {f}")
    with open(f, mode='wt', encoding='utf-8') as myfile:
        for vcf_line in new_vcf:
            myfile.write(vcf_line + '\n')
    logging.debug(f"Sorting, bgzipping {f} -> {out_file}")
    subprocess.run(f"gsort {f} {genome_file} | bgzip -cf > {out_file}", shell=True)
    logging.debug(f"Tabix indexing {out_file}")
    tabix(out_file)
    os.remove(f)

def create_genome_file(vcf_file, genome_file):
    """
    gsort (https://github.com/brentp/gsort) requires a '.genome'
    file to tell it how to sort the vcf records. This function creates a
    '.genome' file in the same order as the contig lines in the vcf header.

    :param vcf_file: a vcf file with the correct contig names
    :param genome_file: the resulting .genome file
    :return: None. this creates a '.genome' file
    """
    
    vcf = vcf.Reader(filename=vcf_file)
    with open(genome_file, mode='wt', encoding='utf-8') as genome_file:
        for contig in vcf.contigs.keys():
            genome_file.write('\t'.join(vcf.contigs[contig]))

def check_dependency(dep, exit=True):
    """
    Check if a dependency exists.

    >>> check_dependency("ls")
    >>> check_dependency("freebayes")
    :param dep: name of the dependency
    :param exit: If True, then if the dependency isn't found, the session will exit.
    :return: True/False if dependency was found.
    """
    found = False
    try:
        res = subprocess.run(['which', dep], capture_output=True, check=True, encoding="UTF8")
        found = True
        logging.info("Found dependency: " + str(res.stdout).strip())
    except subprocess.CalledProcessError:
        logging.error("Missing dependency: " + dep)
        if exit:
            sys.exit(1)
    return found


def check_dependencies(f='verchew.ini'):
    """
    Check that all of the dependencies that are listed in the INI format exist.
    :param f: an INI formatted file. see verchew [https://verchew.readthedocs.io/en/latest/]
    :return:
    """
    config = configparser.ConfigParser()
    config.read(f)
    for section in config.sections():
        dependency = config[section]['cli']
        check_dependency(dependency)