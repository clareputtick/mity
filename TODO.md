# call
* define INFO:tier based on code in report

# normalise
* migrate to pyvcf where possible
* Should document how we combine variants

# merge
* migrate to `pyvcf` where possible (started in `dev/merge2.py`)
* force the same sample names in the same order in mity and HC?
* are there any VCF format 4.1 vs 4.2 fields that clash. eg the type of the variable changes?

# report
* mity report fails if mity normalise hasn't been run, so consider making this mandatory & dropping `mity normalise` 
* doubles up some lines in the report. i've noticed this happen if the same chr-pos-ref-alt is in `phylotree_haplotype` twice
* Check if report works with multiple VCFs
* update mitomap annotations 

# misc
* use logging.info, logging.debug, logging.warning, logging.error where possible
* hg19 support: mity call and merge should be ok with hg19's chrM, but mity report 
uses annotations with GRCh37 coordinates, and thus should fail. How much do we try 
to support hg19 then? GRCh38 and GRCh37 are the same length.

# installation (pre-submission)
* use distutils to create a package, and register this with pip install dev - DONE
* register package with main pip repo - PENDING
* create a docker image - DONE
* update docker image once in main pip repo - PENDING

# GitHub (pre-submission)
* merge branch back to master
* push to KCCG

# DNAnexus
* migrate app code to use the latest mity. either via an asset, or Docker image.
* send app code to GitHub

# Testing (pre-submission)
* seek independent users to test this from scratch
