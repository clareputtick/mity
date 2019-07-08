
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

# MJC
<!-- ## create fasta file
echo -e "MT\t1\t16569" > human_g1k_v37_decoy.MT.bed
samtools faidx $HOME/data/gatk-resource-bundle/2.8/b37/human_g1k_v37_decoy.fasta -o human_g1k_v37_decoy.MT.fa -r human_g1k_v37_decoy.MT.bed -->

## local tests
    mkdir test
    dx download project-F5bzbyQ0Jzv62xBg0P0Q5ygq:/test/kccg-mity-call/inputs/* -o test/
    mity call --reference $B37D5 --prefix test1 test/A1.dedup.realigned.recalibrated.chrMT.bam
    
    freebayes -f /Users/marcow/data/gatk-resource-bundle/2.8/b37/human_g1k_v37_decoy.fasta -b test/A1.dedup.realigned.recalibrated.chrMT.bam -r MT:1-1000 --min-mapping-quality 30 --min-base-quality 20 --min-alternate-fraction 0.5 --min-alternate-count 4 --ploidy 2 --vcf unnormalised.vcf.gz
    
    mity normalise --vcf test1.dedup.realigned.recalibrated.chrMT.mity.vcf.gz --outfile test1.dedup.realigned.recalibrated.chrMT.mity.norm.vcf.gz
