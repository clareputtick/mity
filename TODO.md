# call
* edit the VCF header to include the mity call command
* edit the VCF header to add freebayes_ prefix to any of the important freebayes metadata

# normalise
* migrate to pyvcf where possible

# merge
* migrate to pyvcf where possible (started in dev/merge2.py)
* force the same reference genome in mity and HC?
* force the same sample names in the same order in mity and HC?
* are there any VCF format 4.1 vs 4.2 fields that clash. eg the type of the variable changes?

# report
* L63-L104 is repetitive
* VEP splitting should be in a function
* L195-228 is repetitive
* some more examples in here that could be streamlined and processed over a list of keys
* L823-891 could just be saved in a text file and loaded in as a one-liner
* L903-906 is too repetetive: iterate over an array of fields for int64 vs float64


# misc
* there's no entrypoint to util.create_genome_file, thus gsort may fail
* use logging.info, logging.debug, logging.warning, logging.error where possible

# installation
* use distutils to create a package, and register this with pip install
* create a docker image

# GitHub
* improve documentation
* add INSTALL.md for those that want to install manually
* include DNAnexus app code
* ensure there is example usage

# DNAnexus
* migrate app code to use the latest mity. either via an asset, or Docker image.

# Testing
* seek independent users to test this from scratch
