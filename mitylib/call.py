"""Mitochondrial variant calling."""
import sys
import subprocess
import logging
import os.path
from .util import tabix, check_missing_file, create_prefix, bam_get_mt_contig
from .normalise import do_normalise as vcfnorm

def do_call(bam_files, reference, prefix=None, min_mq=30, min_bq=24,
            min_af=0.01, min_ac=4, p=0.002, normalise=True, 
            out_folder_path=".", region=None):
    """
    Run mity call.
    :param bam_files: a list of bam_files
    :param reference: the path to the reference genome file (fasta format)
    :param prefix: The result filename prefix. If None, then the first bam_file prefix
    will be used
    :param min_mq: minimum mapping quality threshold. default: 30
    :param min_bq: minimum base quality threshold. default: 24
    :param min_af: minimum heteroplasmy, aka the minimum fraction of alt reads vs total reads.
    scale [0,1]; default 0.005
    :param min_ac: minimum number of alternative reads to support a variant. default: 4
    :param p: the noise threshold. default 0.002
    :param normalise:
    :param out_folder_path: the folder to store the results within. default: .
    :param region: Which region to analyse? If None, then the whole MT will be analysed.
    :return:
    """
    bam_files = bam_files[0] 
    #####
    # Checks
    #####

    if len(bam_files) > 1 and prefix is None:
        raise ValueError(
                "If there is more than one bam file, --prefix must be set")
    
    check_missing_file(bam_files, die=True)
    prefix = create_prefix(bam_files[0], prefix)
    
    if not os.path.exists(out_folder_path):
        os.makedirs(out_folder_path)

    output_file_name = os.path.join(out_folder_path, prefix + ".mity.vcf.gz")
    unnormalised_vcf_path = os.path.join(out_folder_path, prefix + ".unnormalised.vcf.gz")

    bam_str = " ".join(['-b ' + bam_file for bam_file in bam_files])

    if region is None:
        region = bam_get_mt_contig(bam_files[0], as_string=True)

    # this is the mity command that will be printed in the VCF header
    # I'm not sure how to get the exact line from the terminal, so I remake it here
    mity_cmd = '##commandline="mity call --reference ' + str(reference) + ' --prefix ' + prefix + ' --min-mapping-quality ' + str(min_mq) + ' --min-base-quality ' + str(min_bq) + ' --min-alternate-fraction ' + str(min_af) + ' --min-alternate-count ' + str(min_ac)  + ' --out-folder-path ' + str(out_folder_path)

    if region is not None:
        mity_cmd = mity_cmd + ' --region ' + region

    if normalise:
        mity_cmd = mity_cmd + ' --normalise --p ' + str(p)

    mity_cmd = mity_cmd + '"'
    # print(mity_cmd)

    mity_cmd = mity_cmd.replace("/", "\/")

    # print(mity_cmd)

    # replace the freebayes command line with the mity command
    sed_cmd = "sed 's/^##commandline=.*/" + mity_cmd + "/'"
    # print(sed_cmd)

    freebayes_call = ('freebayes -f {} {} '
                      '--min-mapping-quality {} '
                      '--min-base-quality {} '
                      '--min-alternate-fraction {} '
                      '--min-alternate-count {} '
                      '--ploidy 2 '
                      '--region {} '
                      ).format(reference, bam_str, min_mq, min_bq, min_af, min_ac, region)

    freebayes_call = freebayes_call + ('| {} | bgzip > {} ').format(sed_cmd, unnormalised_vcf_path)

    logging.info("Running FreeBayes in sensitive mode")
    # print(freebayes_call)
    
    subprocess.run(freebayes_call, shell=True)
    if os.path.isfile(unnormalised_vcf_path):
        logging.debug("Finished running FreeBayes")
    
    if normalise:
        logging.debug("Normalising and Filtering variants")
        try:
            vcfnorm(vcf=unnormalised_vcf_path, out_file=output_file_name, p=p)
        finally:
            os.remove(unnormalised_vcf_path)
    else:
        os.rename(unnormalised_vcf_path, output_file_name)
        tabix(output_file_name)
