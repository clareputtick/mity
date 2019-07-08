#!/Users/clareputtick/anaconda/bin/python3

import sys
import argparse
import subprocess

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run mity-call')
    parser.add_argument('bam', action = 'append', nargs='+', help = 'BAM files to run the analysis on.')
    parser.add_argument('--reference', action='store', required=True)
    parser.add_argument('--prefix', action='store', help = 'Output files will be named with PREFIX')
    parser.add_argument('--sites', action='store', help = 'Force freebayes to call variants at these sites. Should be vcf.gz.')
    parser.add_argument('--min-mapping-quality', action='store', type = int, default = 30, help = 'Exclude alignments from analysis if they have a mapping quality less than MIN_MAPPING_QUALITY. Default: 30')
    parser.add_argument('--min-base-quality', action='store', type = int, default = 20, help = 'Exclude alleles from analysis if their supporting base quality is less than MIN_BASE_QUALITY. Default: 20')
    parser.add_argument('--min-alternate-fraction', action='store', type = float, default = 0.5, help = 'Require at least MIN_ALTERNATE_FRACTION observations supporting an alternate allele within a single individual in the in order to evaluate the position. Default: 0.0001, range = [0,1]')
    parser.add_argument('--min-alternate-count', action='store', type = int, default = 4, help = 'Require at least MIN_ALTERNATE_COUNT observations supporting an alternate allele within a single individual in order to evaluate the position. Default: 4')
    parser.add_argument('--ploidy', action='store', type = int, default = 2, help = 'Expected ploidy of the sample. Default: 2.')
    # parser.add_argument('--parallel', action='store_true', help = 'Run freebayes in parallel.')
    # parser.add_argument('--ncpu', action='store', type = int, help = 'Number of CPUs to use when running in parallel.')

    
    args = parser.parse_args()
    # the only positional arguments are bam files
    bam_files = args.bam[0]
    # print(args.bam)
    # sys.exit()

    #####
    # Checks
    #####
    # TODO: check bam and bai match, or just input bam and find bai in directory
    # TODO: check sites parameters only has 8 columns
    # TODO: need to check what the mitochondria is called in the bam header? 

    # Check if no prefix set that there is only one bam
    if args.prefix is None and len(bam_files) > 1:
      # print(bam_files)
      # print(len(bam_files))
      sys.exit("Error: If more than one bam, --prefix must be set")

    if args.prefix is None:
      # then only one bam file
      prefix=bam_files[0].split(".")[0]
    else:
      prefix=args.prefix
    
    # make outfile string
    file_string = []
    for b in bam_files:
      f_string = ".".join(b.split('.')[1:-1])
      f_string=f_string + '.mity.vcf.gz'
      file_string.append(f_string)

    # check all the bam files have the same file_string
    if len(set(file_string)) > 1:
      sys.exit("Don't know what to do here")
    file_string = file_string[0]

    output_file_name = prefix + "." + file_string
    
    freebayes_call="freebayes -f " + args.reference
    for b in bam_files:
      freebayes_call = freebayes_call + ' -b ' + b
    
    freebayes_call = freebayes_call + " --min-mapping-quality " + str(args.min_mapping_quality)
    freebayes_call = freebayes_call + " --min-base-quality " + str(args.min_base_quality)
    freebayes_call = freebayes_call + " --min-alternate-fraction " + str(args.min_alternate_fraction)
    freebayes_call = freebayes_call + " --min-alternate-count " + str(args.min_alternate_count)
    freebayes_call = freebayes_call + " --ploidy " + str(args.ploidy)

    freebayes_call = freebayes_call + " | bgzip > unnormalised.vcf.gz"

    # run freebayes
    # https://stackoverflow.com/questions/89228/calling-an-external-command-in-python
    # subprocess.run('freebayes -f ../../../hs37d5.fasta-index/genome.fa -b A7.dedup.realigned.recalibrated.chrMT.bam | bgzip > test.vcf.gz', shell=True )
    sys.stderr.write("Running FreeBayes\n")
    # print(freebayes_call)
    subprocess.run(freebayes_call, shell=True )
    
    # run mity normalise
    mity_normalise_call = "/Users/clareputtick/Garvan/mt_apps/kccg-mity-call/resources/home/dnanexus/mity_normalise.py --vcf unnormalised.vcf.gz |  sort -k1,1 -k2,2n | bgzip > " + output_file_name
    sys.stderr.write("Normalising variants\n")
    # print(mity_normalise_call)
    subprocess.run(mity_normalise_call, shell=True )
    
    tabix_call = "tabix " + output_file_name
    # print(tabix_call)
    subprocess.run(tabix_call, shell=True )
    subprocess.run('rm unnormalised.vcf.gz', shell=True )
