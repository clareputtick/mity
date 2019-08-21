# mity
mity is a bioinformatic analysis pipeline designed to call mitochondrial SNV and INDEL variants from Whole Genome Sequencing (WGS) data. mity can:
* identify very low-heteroplasmy variants, even <1% heteroplasmy when there is sufficient read-depth (eg >1000x)
* filter out common artefacts that arise from high-depth sequencing
* easily integrate with existing nuclear DNA analysis pipelines (mity merge)
* provide an annotated report, designed for clinicians and researchers to interrogate


# Usage
    mity -h

# Dependencies
* freebayes >= 1.2.0
* bgzip + tabix
* gsort (https://github.com/brentp/gsort)
* python3 (tested on 3.7.4)
* pyvcf
* xlsxwriter
* pandas

# Example Usage
This is an example of calling variants in the Ashkenazim Trio.

First make sure mity is in your PATH variable.

```bash
PATH="PATH_TO_MITY_FOLDER:${PATH}"
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
test_in/HG002.hs37d5.2x250.small.MT.RG.bam \
test_in/HG003.hs37d5.2x250.small.MT.RG.bam \
test_in/HG004.hs37d5.2x250.small.MT.RG.bam 
```

This should create test_out/ashkenazim.mity.vcf.gz and test_out/ashkenazim.mity.vcf.gz.tbi


# Acknowledgements
We thank the Kinghorn Centre for Clinical Genomics and collaborators, who helped
with feedback for running mity. 
We thank Eric Talevich who's CNVkit helped us structure mity as a package
