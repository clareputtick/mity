"""Mitochondrial variant calling."""
import sys
import subprocess
import logging
import os.path
from .util import tabix
from .util import check_missing_file
from .normalise import do_normalise as vcfnorm

# from . import normalise ## @TODO

def do_call(bam_files, reference, prefix=None, min_mq=30, min_bq=20,
            min_af=0.5, min_ac=4, ploidy=2, normalise=True):
    bam_files = bam_files[0]  ## @TODO check this still works with >1 BAM file
    
    #####
    # Checks
    #####
    # TODO: need to check what the mitochondria is called in the bam header? 
    

    if len(bam_files) > 1 and prefix is None:
        raise ValueError(
                "If there is more than one bam file, --prefix must be set")
    
    check_missing_file(bam_files, die=True)
    
    prefix = [bam_files[0].split(".")[0], prefix][prefix is not None]
    
    # make outfile string
    file_string = []
    for b in bam_files:
        f_string = ".".join(b.split('.')[1:-1]) + '.mity.vcf.gz'
        file_string.append(f_string)
    
    # check all the bam files have the same file_string
    if len(set(file_string)) > 1:
        sys.exit("Don't know what to do here")
    file_string = file_string[0]
    
    output_file_name = prefix + "." + file_string
    
    bam_str = " ".join(['-b ' + bam_file for bam_file in bam_files])
    
    region = "MT:1-16569"  # @TODO parse chrom name & length from the BAM header
    region = "MT:1-500"  # @TODO delete this debugging sub-region analysis
    freebayes_call = (f'freebayes -f {reference} {bam_str} -r {region} ' +
                      f'--min-mapping-quality {min_mq} ' +
                      f'--min-base-quality {min_bq} ' +
                      f'--min-alternate-fraction {min_af} ' +
                      f'--min-alternate-count {min_ac} ' +
                      f'--ploidy {ploidy} | bgzip > unnormalised.vcf.gz'
                      )
    # run FreeBayes
    logging.info("Running FreeBayes in sensitive mode")
    print(freebayes_call)
    subprocess.run(freebayes_call, shell=True)
    if os.path.isfile('unnormalised.vcf.gz'):
        logging.info("Finished running FreeBayes")
    
    if normalise:
        logging.info("Normalising and FILTERing variants")
        try:
            vcfnorm('unnormalised.vcf.gz', output_file_name)
        finally:
            os.remove('unnormalised.vcf.gz')
    else:
        os.rename("unnormalised.vcf.gz", output_file_name)
        tabix(output_file_name)


