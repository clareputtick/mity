
# Run version 1

## Call variants in the mitochondria
python3 call.py --reference genome.MT.fa --min-alternate-fraction 0.5 --prefix 15F00007 SYD-40100730.dedup.realigned.recalibrated.chrMT.bam SYD-40100741.dedup.realigned.recalibrated.chrMT.bam

## merge mity vcf and hc vcf
python3 merge.py --mity 15F00007.dedup.realigned.recalibrated.chrMT.mity.vcf.gz --hc 15F00007.hc.small.vcf.gz

## make mity report
python3 report.py 15F00007.dedup.realigned.recalibrated.chrMT.mity.vcf.gz

###
# Run version 2
##

## Call variants in the mitochondria
python3 call.py --reference genome.MT.fa --min-alternate-fraction 0.5  SYD-40100741.dedup.realigned.recalibrated.chrMT.bam

python3 call.py --reference genome.MT.fa --min-alternate-fraction 0.5 SYD-40100730.dedup.realigned.recalibrated.chrMT.bam 

## make mity report
python3 report.py SYD-40100730.dedup.realigned.recalibrated.chrMT.mity.vcf.gz SYD-40100741.dedup.realigned.recalibrated.chrMT.mity.vcf.gz --prefix 15F00007