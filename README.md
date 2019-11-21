# mity
`mity` is a bioinformatic analysis pipeline designed to call mitochondrial SNV and INDEL variants from Whole Genome Sequencing (WGS) data. `mity` can:
* identify very low-heteroplasmy variants, even <1% heteroplasmy when there is sufficient read-depth (eg >1000x)
* filter out common artefacts that arise from high-depth sequencing
* easily integrate with existing nuclear DNA analysis pipelines (mity merge)
* provide an annotated report, designed for clinicians and researchers to interrogate

# Usage
    mity -h

# Dependencies
* python3 (tested on 3.7.4)
* freebayes >= 1.2.0
* bgzip + tabix
* gsort (https://github.com/brentp/gsort)
* pyvcf
* xlsxwriter
* pandas

# Installation
Installation instructions via Docker, pip, or manually are available in [INSTALL.md](https://github.com/KCCG/mity/blob/master/INSTALL.md)

# Example Usage
This is an example of calling variants in the Ashkenazim Trio.

## mity call
First run `mity call` on three MT BAMs provided in [mity/test_in](https://github.com/KCCG/mity/blob/master/test_in).

We recommend always using `--normalise`, or `mity report` won't work:
```bash
mity call \
--prefix ashkenazim \
--out-folder-path test_out \
--region MT:1-500 \
--normalise \
test_in/HG002.hs37d5.2x250.small.MT.RG.bam \
test_in/HG003.hs37d5.2x250.small.MT.RG.bam \
test_in/HG004.hs37d5.2x250.small.MT.RG.bam 
```
This will create `test_out/normalised/ashkenazim.mity.vcf.gz` (and tbi file).

## mity report

We can create a `mity report` on the normalised VCF:
```bash
mity report \
--prefix ashkenazim \
--min_vaf 0.01 \
--out-folder-path test_out \
test_out/ashkenazim.mity.vcf.gz
```
This will create: `test_out/ashkenazim.annotated_variants.csv` and `test_out/ashkenazim.annotated_variants.xlsx`.

## mity normalise
High-depth sequencing and sensitive variant calling can create many variants with more than 2 alleles, and in some
cases, joins two nearby variants separated by shared `REF` sequence into a multi-nucleotide polymorphism 
as discussed in the manuscript. Here, variant normalisation relates to decomposing the multi-allelic variants and 
where possible, splitting multi-nucleotide polymorphisms into their cognate smaller variants. At the time of writing,
all variant decomposition tools we used failed to propagate the metadata in a multi-allelic variant to the split
variants which caused problems when reporting the quality scores associated with each variant.
  
Technically you can run `mity call` and `mity normalise` separately, but since `mity report` requires a normalised 
vcf file, we recommend running `mity call --normalise`. 

## mity merge
You can merge a nuclear vcf.gz file and a mity.vcf.gz file thereby replacing the MT calls from the nuclear VCF (
presumably from a caller like HaplotypeCaller which is not able to sensitively call mitochondrial variants) with
the calls from `mity`.

```bash
mity merge \
--prefix ashkenazim \
--mity_vcf test_out/ashkenazim.mity.vcf.gz \
--nuclear_vcf todo-create-example-nuclear.vcf.gz
```

# Recommendations for interpreting the report
Assuming that you are looking for a pathogenic variant underlying a patient with a rare genetic disorder potentially 
caused by a Mitochondrial mutation, then we recommend the following strategy:
1. tier 1 or 2 variants included in the 'commercial_panels' column 
2. tier 1 or 2 variants that match the clinical presentation and the phenotype in 'disease_mitomap', preferably 
those that are annotated with Confirmed evidence in the 'status_mitomap' column
3. exclude common variants: anything linked to 'phylotree_haplotype', high 'phylotree_haplotype', high 
'MGRB_frequency', high 'GenBank_frequency_mitomap'.
4. consider any remaining tier 1 or 2 variants that may have a predicted impact on tRNA
5. consider any remaining variants with high numbers of 'variant_references_mitomap'
5. if you have analysed multiple family members, consider variants who's level of 'variant_heteroplasmy' match the
disease burden 

# Acknowledgements
We would like to thank:
* The Kinghorn Centre for Clinical Genomics and collaborators, who helped with feedback for running `mity`.
* The Genome in a Bottle consortium for providing the test data used here 
* Eric Talevich who's CNVkit helped us structure `mity` as a package
* Erik Garrison for developing `FreeBayes` and his early feedback in optimising `FreeBayes` for sensitive variant detection.
* Brent Pederson for developing `gsort`