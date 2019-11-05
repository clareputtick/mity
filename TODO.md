# call
* Check that hg19 reference works

# normalise
* migrate to pyvcf where possible
* Should document how we combine variants

# merge
* migrate to `pyvcf` where possible (started in `dev/merge2.py`)
* force the same sample names in the same order in mity and HC?
* are there any VCF format 4.1 vs 4.2 fields that clash. eg the type of the variable changes?

# report
* mity report fails if mity normalise hasn't been run, so consider making this mandatory & dropping `mity normalise` 
* CRITICAL: currently broken:
    File "/usr/local/lib/python3.7/site-packages/mitylib/report.py", line 353, in split_header_variants
        col_names = header[-1]
    IndexError: list index out of range
>> This is because the mity-report code assumes certain fields are in the VCF, which aren't there if mity-normalise hasnt been run. I think the best way to fix this is to rewrite the code so that it doesn't assume fields in the VCF. This could then be a more useful script as it would work with all VCFs. But I'm going to leave it for now, and we can say that mity-report only works with normalised VCFs.

* L63-L104 is repetitive
* VEP splitting should be in a function
* L195-228 is repetitive
* some more examples in here that could be streamlined and processed over a list of keys
* L823-891 could just be saved in a text file and loaded in as a one-liner
* L903-906 is too repetitive: iterate over an array of fields for int64 vs float64
* Check if report works with multiple VCFs
* need to add output folder - currently saves in current directory - Done

# misc
* use logging.info, logging.debug, logging.warning, logging.error where possible
* hg19 support: mity call and merge should be ok with hg19's chrM, but mity report 
uses annotations with GRCh37 coordinates, and thus should fail. How much do we try 
to support hg19 then? GRCh38 and GRCh37 are the same length.

# installation (pre-submission)
* use distutils to create a package, and register this with pip install - to finalise
* create a docker image - DONE

# GitHub (pre-submission)
* CRITICAL: improve documentation
* CRITICAL: update INSTALL.md
* CRITICAL: ensure there is example usage

# DNAnexus
* migrate app code to use the latest mity. either via an asset, or Docker image.
* send app code to GitHub

# Testing (pre-submission)
* seek independent users to test this from scratch
