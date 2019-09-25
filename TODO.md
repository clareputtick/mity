# call
* CRITICAL: makes an odd, very long filename - DONE
* CRITICAL: edit the VCF header to include the mity call command
* CRITICAL: edit the VCF header to add freebayes_ prefix to any of the important freebayes metadata
* CRITICAL: If the freebayes command (or any other subprocess.run commands) fail, mity should fail
* At the moment we use gsort or bcftools tools to sort the normalised vcf. Given that it is always chromosome MT, this should be able to be done with python in the normalise script, meaning one less tool to download.
* Check that the different references work
* Freebayes assumes that the BAMs have a RG, which is where it gets the sample name from. If there is no @RG line, freebayes just outputs a single "unknown" sample, even if you input more than one sample. We should add a check to mity, to check that the BAM header has a read group line. This could just check that there is a line in the BAM header starting with @RG.
* I dont think the error for if the bams dont exist is working:
```bash
mity call --prefix ashkenazim --out-folder-path test_out --min-alternate-fraction 0.5 --normalise bam_that_doesnt_exist.bam
```

# normalise
* mity report fails if mity normalise hasn't been run, so consider making this mandatory & dropping `mity normalise` 
* migrate to pyvcf where possible
* CRITICAL: i've added the 'p' noise floor threshold, but can't see where the new QUAL is calculated. - Done
* CRITICAL: When there is more than one sample I dont think the sample names are coming out properly - Done
* Make the test bams contain lines that have "repeated positions". To do this fun fb on entire bams and see which regions would give a repeated position. - Done
* Should probably detail how we combine variants
* If the number of reads supporting the variant and the depth are the same, we get a divide by zero error, because the binomial cdf gives 1, and we end up with log10(1-1). We need to decide what q should be in this case - Inf or a large number?

# merge
* migrate to pyvcf where possible (started in dev/merge2.py)
* force the same sample names in the same order in mity and HC?
* are there any VCF format 4.1 vs 4.2 fields that clash. eg the type of the variable changes?

# report
* CRITICAL: currently broken:
    File "/usr/local/lib/python3.7/site-packages/mitylib/report.py", line 353, in split_header_variants
        col_names = header[-1]
    IndexError: list index out of range

* L63-L104 is repetitive
* VEP splitting should be in a function
* L195-228 is repetitive
* some more examples in here that could be streamlined and processed over a list of keys
* L823-891 could just be saved in a text file and loaded in as a one-liner
* L903-906 is too repetitive: iterate over an array of fields for int64 vs float64

# misc
* use logging.info, logging.debug, logging.warning, logging.error where possible
* hg19 support: mity call and merge should be ok with hg19's chrM, but mity report 
uses annotations with GRCh37 coordinates, and thus should fail. How much do we try 
to support hg19 then? GRCh38 and GRCh37 are the same length.

# installation (pre-submission)
* use distutils to create a package, and register this with pip install
* create a docker image

# GitHub (pre-submission)
* CRITICAL: improve documentation
* CRITICAL: update INSTALL.md
* CRITICAL: ensure there is example usage

# DNAnexus
* migrate app code to use the latest mity. either via an asset, or Docker image.
* send app code to GitHub

# Testing (pre-submission)
* seek independent users to test this from scratch
