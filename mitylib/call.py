"""Mitochondrial variant calling."""
import sys
import subprocess
import logging
import os.path


# from . import normalise ## @TODO

def do_call(bam_files, reference, prefix=None, min_mq=30, min_bq=20,
            min_af=0.5, min_ac=4, ploidy=2):
    bam_files = bam_files[0]  ## @TODO check this still works with >1 BAM file
    
    #####
    # Checks
    #####
    # TODO: check bam and bai match, or just input bam and find bai in directory
    # TODO: check sites parameters only has 8 columns
    # TODO: need to check what the mitochondria is called in the bam header? 
    
    # Check if no prefix set that there is only one bam
    if len(bam_files) == 0:
        raise ValueError("At least one BAM file must be supplied (-b)")
    if len(bam_files) > 1 and prefix is None:
        raise ValueError(
                "If there is more than one bam file, --prefix must be set")
    
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
    
    bam_str = ['-b ' + bam_file for bam_file in bam_files]
    
    region = "MT:1-16569"  # @TODO parse chrom name & length from the BAM header
    region = "MT:1-500"  # @TODO delete this debugging sub-region analysis
    freebayes_call = (f'freebayes -f {reference} {bam_str} -r {region} ' +
                      f'--min-mapping-quality {min_mq} ' +
                      f'--min-base-quality {min_bq} ' +
                      f'--min-alternate-fraction {min_af} ' +
                      f'--min-alternate-count {min_ac} ' +
                      f'--ploidy {ploidy} --vcf unnormalised.vcf.gz'
                      )
    # run FreeBayes
    logging.info("Running FreeBayes in sensitive mode")
    print(freebayes_call)
    subprocess.run(freebayes_call, shell=True)
    if os.path.isfile('unnormalised.vcf.gz'):
        logging.info("Finished running FreeBayes")
    
    # run mity normalise
    # norm_vcf = normalise.do_call('unnormalised.vcf.gz', 'MT', 
    # output_file_name)
    mity_normalise_call = "~/src/DNANexus_Applets/kccg-mity-call/resources" \
                          "/home/dnanexus/mity_normalise.py --vcf " \
                          "unnormalised.vcf.gz | sort -k1,1 -k2,2n | bgzip > " \
                          "" + output_file_name
    sys.stderr.write("Normalising variants\n")
    # print(mity_normalise_call)
    subprocess.run(mity_normalise_call, shell=True)
    
    tabix_call = "tabix " + output_file_name
    # print(tabix_call)
    subprocess.run(tabix_call, shell=True)
    subprocess.run('rm unnormalised.vcf.gz', shell=True)
