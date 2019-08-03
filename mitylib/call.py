"""Mitochondrial variant calling."""
import sys
import subprocess
import logging
import os.path
from .util import tabix
from .util import check_missing_file
from .util import create_prefix
from .normalise import do_normalise as vcfnorm

def do_call(bam_files, reference, prefix=None, min_mq=30, min_bq=24,
            min_vaf=0.005, min_ac=4, p=0.002, normalise=True, out_folder_path="/Users/putticc/Projects/mity/test_out"):
    """
    Run mity call.
    :param bam_files: a list of bam_files
    :param reference: the path to the reference genome file (fasta format)
    :param prefix: The result filename prefix. If None, then the first bam_file prefix
    will be used
    :param min_mq: minimum mapping quality threshold. default: 30
    :param min_bq: minimum base quality threshold. default: 24
    :param min_vaf: minimum heteroplasmy, aka the minimum fraction of alt reads vs total reads.
    scale [0,1]; default 0.005
    :param min_ac: minimum number of alternative reads to support a variant. default: 4
    :param p: the noise threshold. default 0.002
    :param normalise:
    :return:
    """
    bam_files = bam_files[0] 
    #####
    # Checks
    #####
    # TODO: need to check what the mitochondria is called in the bam header? 
    

    if len(bam_files) > 1 and prefix is None:
        raise ValueError(
                "If there is more than one bam file, --prefix must be set")
    
    check_missing_file(bam_files, die=True)
    prefix = create_prefix(bam_files[0], prefix)
    
    if prefix is None:
        prefix=os.path.basename(bam_files[0]).replace(".bam", "") # @TODO need to check that the bam ends in .bam?


    output_file_name = out_folder_path + "/" + prefix + ".mity.vcf.gz"
    
    bam_str = " ".join(['-b ' + bam_file for bam_file in bam_files])
    unnormalised_vcf_path = out_folder_path + "/" + prefix + ".unnormalised.vcf.gz"

    region = "MT:1-16569"  # @TODO parse chrom name & length from the BAM header
    region = "MT:1-500"  # @TODO delete this debugging sub-region analysis
    freebayes_call = ('freebayes -f {} {} -r {} '
                      '--min-mapping-quality {} '
                      '--min-base-quality {} '
                      '--min-alternate-fraction {} '
                      '--min-alternate-count {} '
                      '--ploidy 2 | bgzip > {} '
                      ).format(reference, bam_str, region, min_mq, min_bq, min_vaf, min_ac, unnormalised_vcf_path)

    logging.info("Running FreeBayes in sensitive mode")
    print(freebayes_call)
    subprocess.run(freebayes_call, shell=True)
    if os.path.isfile(unnormalised_vcf_path):
        logging.info("Finished running FreeBayes")
    
    if normalise:
        logging.info("Normalising and FILTERing variants")
        try:
            vcfnorm(unnormalised_vcf_path, output_file_name, p)
        finally:
            os.remove(unnormalised_vcf_path)
    else:
        os.rename(unnormalised_vcf_path, output_file_name)
        tabix(output_file_name)


