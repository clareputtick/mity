# mity
@TODO

# Usage
    mity -h

# Dependencies
* freebayes >= 1.2.0
* bgzip + tabix
* gsort (https://github.com/brentp/gsort)
* python3 (tested on 3.7.3)
* pyvcf
* xlsxwriter
* pandas

# Example Usage
This is an example of calling variants in the Ashkenazim Trio.

First make sure mity is in your PATH variable.

```bash
PATH="PATH_TO_MITY:${PATH}"
export PATH
```

## mity-call
First run mity-call on three MT BAMs provided in mity/test_in

```bash
mity call \
--prefix ashkenazim \
--out-folder-path test_out \
--min-alternate-fraction 0.5 \
--region MT:1-500 \
--normalise \
test_in/HG002.hs37d5.2x250.MT.bam \
test_in/HG003.hs37d5.2x250.MT.bam \
test_in/HG004.hs37d5.2x250.MT.bam 
```

This should create...

@TODO: 
getting error tabix -f test_out/ashkenazim.mity.vcf.gz
[E::hts_idx_push] Unsorted positions on sequence #82: 16569 followed by 263

# Acknowledgements
We thank the Kinghorn Centre for Clinical Genomics and collaborators, who helped
with feedback for running mity. 
We thank Eric Talevich who's CNVkit helped us structure mity as a package
