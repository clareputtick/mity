# call
* define INFO:tier based on code in report
* add CRAM support

# normalise
* migrate to `pyvcf` where possible
* Should document how we combine variants

# merge
* migrate to `pyvcf` where possible (started in `dev/merge2.py`)
* force the same sample names in the same order in mity and HC?
* are there any VCF format 4.1 vs 4.2 fields that clash. eg the type of the variable changes?

# report
* mity report fails if mity normalise hasn't been run, so consider making this mandatory & dropping `mity normalise` 
* Check if report works with multiple VCFs
* update mitomap annotations 

# misc
* use logging.info, logging.debug, logging.warning, logging.error where possible
* hg19 support: after running BLAST, there are more sequence differences than initially appreciated.
 This is scheduled for deprecation

# DNAnexus
* migrate app code to use the latest mity. either via an asset, or Docker image.
