Cleanup the mitmap annotation file

* rename some column headers:
    perl -pi -e 's/variant_references_mitomap/num_references_mitomap/' mitomap_panel_annotations.csv
    perl -pi -e 's/R\._mitomap/RNA_mitomap/g' mitomap_panel_annotations.csv
    perl -pi -e 's/GB_frequency_mitomap/GenBank_frequency_mitomap/g' mitomap_panel_annotations.csv
* replace the NAs
    perl -pi~ -e 's/NA/./g' mitomap_panel_annotations.csv
* 7471insC == C7471CC, so merge the annotations from these 2 records
* 7472insC == C7471CC, so delete this
* still need to fix the del POS
* replace the Cfrm
    perl -pi -e 's/Cfrm/Confirmed/g' mitylib/annot/mitomap_panel_annotations.csv

# fix the del records
library(mjcbase)
orig <- read.csv("mitylib/annot/mitomap_panel_annotations.csv", stringsAsFactors=FALSE)
head(orig)
orig$commercial_panels <- ifelse(orig$baylor_panel=="TRUE", "Baylor;", "")
orig$commercial_panels <- paste(orig$commercial_panels, ifelse(orig$common_22_panel=="TRUE", "Common22;", ""), sep="")
orig$commercial_panels <- paste(orig$commercial_panels, ifelse(orig$common_58_panel=="TRUE", "Common58;", ""), sep="")
orig$commercial_panels <- sub(";$", "", orig$commercial_panels)
orig$commercial_panels[orig$commercial_panels==""] <- "."
table(orig$commercial_panels)

orig[orig==".;."] <- "."
orig[orig==".;.;."] <- "."
orig[orig==".;.;.;,"] <- "."
orig[orig==".;.;.;."] <- "."
orig[orig==".;.;.;.;."] <- "."
orig[orig==".;.;.;.;.;."] <- "."
orig[orig=="rR."] <- "rRNA"
orig[orig=="tR."] <- "tRNA"

#write.csv(orig, "mitomap_panel_annotations_fixed.csv")

# fix the del records

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install()
# BiocManager::install('Biostrings')
library(Biostrings)

?fasta.index
MT <- readDNAStringSet("mitylib/reference/hs37d5.MT.fa")

orig$POSminus1 <- sapply(orig$POS-1, function(x) as.character(DNAStringSet(MT, x, width=1)))
idx <- grepl("del", orig$ALT)
sum(idx)
# [1] 405

orig$NEWREF <- ""
orig$NEWREF[idx] <- paste(orig$POSminus1[idx], orig$REF[idx], sep="")
orig$NEWALT[idx] <- orig[idx,"POSminus1"]
head(orig[idx,])
orig$REF[idx] <- orig$NEWREF[idx]
orig$ALT[idx] <- orig$NEWALT[idx]
orig$POS[idx] <- orig$POS[idx] - 1
# test ref using vt normalize, as it will check REF base for us
# #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG002   HG003   HG004

orig <- orig[, !(colnames(orig) %in% c("POSminus1", "NEWREF", "NEWALT"))]
write.csv(orig, "mitomap_panel_annotations_fixed.csv")

vcf <- with(orig,
    data.frame(
        CHROM="MT", 
        POS=POS,
        ID=1:nrow(orig),
        REF=REF,
        ALT=ALT,
        QUAL=".",
        FILTER=".",
        INFO=".")
    )
write.delim(vcf, "mitomap_panel_annotations_fixed.vcf")
# insert:
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=20191107
##reference=/Users/mcowley/src/mity/mitylib/reference/hs37d5.MT.fa							
##contig=<ID=MT,length=16569>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
$ vt normalize -n -r mitylib/reference/hs37d5.MT.fa -o mitomap_panel_annotations_fixed.left.vcf mitomap_panel_annotations_fixed.vcf 2>&1  | grep -c skipped
15
...
reference bases not consistent: MT:513-515  CCA(REF) vs CAC(FASTA)
reference bases not consistent: MT:514-516  AAC(REF) vs ACA(FASTA)
reference bases not consistent: MT:2351-2353  TTA(REF) vs TAA(FASTA)
reference bases not consistent: MT:3105-3106  CC(REF) vs CN(FASTA)
reference bases not consistent: MT:5133-5135  AAA(REF) vs ACT(FASTA)
reference bases not consistent: MT:8275-8283  CCACCCCCT(REF) vs CTCTACCCC(FASTA)
reference bases not consistent: MT:8278-8287  TCCCCCTCTA(REF) vs TACCCCCTCT(FASTA)
reference bases not consistent: MT:8287-8296  TCCCCCTCTA(REF) vs TAGAGCCCAC(FASTA)
reference bases not consistent: MT:9477-9492  TTTTTTCTTCGCAGGA(REF) vs TTTTTTTCTTCGCAGG(FASTA)
reference bases not consistent: MT:11620-11622  TTA(REF) vs TAC(FASTA)
reference bases not consistent: MT:15257-15258  GA(REF) vs AC(FASTA)

# fix this record
reference bases not consistent: MT:16045-16047  T>G(REF) vs TGG(FASTA)
perl -pi -e 's/"T>G",./"T","G"/' mitylib/annot/mitomap_panel_annotations.csv


# Fix these 3 manually
reference bases not consistent: MT:15497-15510  24 bp deletion(REF) vs GCGACCCAGACAAT(FASTA)
reference bases not consistent: MT:15648-15661  18 bp deletion(REF) vs AGCAATAATCCCCA(FASTA)
reference bases not consistent: MT:8927-8937  ACACCT8928d(REF) vs TACACCCCTTA(FASTA)

# Fix 2 records found by vt b
* replace the "24 bp deletion",.
aka 251GDPDNYTL-del258	
grep "24 bp deletion" mitylib/annot/mitomap_panel_annotations.csv | perl -p -e 's/15498,"24 bp deletion",./15497,"GCGACCCAGACAATTATACCCTAG","G"/'
perl -pi -e 's/15498,"24 bp deletion",./15497,"GCGACCCAGACAATTATACCCTAG","G"/' mitylib/annot/mitomap_panel_annotations.csv

* replace the "18 bp deletion"
    grep "18 bp deletion" mitylib/annot/mitomap_panel_annotations.csv 
    15649,"18 bp deletion",.,"MT-CYB",.,.,.,.,"Multisystem Disorder, EXIT","1",.,"-","+","Reported","ILAMIP-del","0",.,.,.,"MT"

    faidx mitylib/reference/hs37d5.MT.fa MT:15648-15667
    >MT:15648-15667
    TAGCAATAATCCCCATCCTC
    
    perl -pi -e 's/15649,"18 bp deletion",./15648,"TAGCAATAATCCCCATCCTC","T"/' mitylib/annot/mitomap_panel_annotations.csv

* ACACCT8928d
Assume this is ACACCT8928del
    grep ACACCT8928d  mitylib/annot/mitomap_panel_annotations.csv
    8928,"ACACCT8928d",.
    hunt...
    >MT:8921-8935
    GCACACCTACACCCC
      ******
     ^m.8922 
     MT 28922 CACACCT C
     
     perl -pi -e 's/8928,"ACACCT8928d",./8922,"CACACCT","C"/' mitylib/annot/mitomap_panel_annotations.csv


* reference bases not consistent: MT:9477-9492  TTTTTTCTTCGCAGGA(REF) vs TTTTTTTCTTCGCAGG(FASTA)
= Myoglobinuria","9480del15","TTTTTCTTCGCAGGA del","FFFAG del"," "," "," …
    faidx mitylib/reference/hs37d5.MT.fa MT:9479-9493
    >MT:9471-9499
    TCAGAAGTTTTTTTCTTCGCAGGATTTTT
             TTTTTCTTCGCAGGAdel
            ^ m.9480

perl -pi -e 's/9479,"TTTTTCTTCGCAGGA","del",/9480,"TTTTTTCTTCGCAGGA","T",/' mitylib/annot/mitomap_panel_annotations.csv

* reference bases not consistent: MT:8275-8283  CCACCCCCT(REF) vs CTCTACCCC(FASTA)
= 8275	CCTCTA	del
    >MT:8271-8290
    ACCCCCTCTACCCCCTCTAG
           CCACCCCCT
          ^m.8278
    perl -pi -e 's/8277,"CACCCCCT","del"/8278,"TCTACCCCCT","T"/' mitylib/annot/mitomap_panel_annotations.csv

* reference bases not consistent: MT:8278-8287  TCCCCCTCTA(REF) vs TACCCCCTCT(FASTA)
= 8280	CCCCCTCTA	del == CCCCCTCTA8281d
    >MT:8271-8290
    ACCCCCTCTACCCCCTCTAG
     CCCCCTCTA
    ^m.8271
    perl -pi -e 's/8272,"CCCCCTCTA","del"/8271,"ACCCCCTCTA","A"/'  mitylib/annot/mitomap_panel_annotations.csv
    
* reference bases not consistent: MT:8287-8296  TCCCCCTCTA(REF) vs TAGAGCCCAC(FASTA)
    faidx mitylib/reference/hs37d5.MT.fa MT:8281-8299
    >MT:8281-8299
    CCCCCTCTAGAGCCCACTG
    CCCCCTCTA
   ^m.8280
   
   # the 8281 will get adjusted to be correct, but will clash with the 8280 record that i'm about to fix
   $ grep CCCCCTCTA mitylib/annot/mitomap_panel_annotations.csv
   8280,"CCCCCTCTA","del","MT-NC7","84","non-coding","-","-",.,.,.,.,.,.,.,"0",.,.,.,"MT"
   8281,"CCCCCTCTA","del","MT-NC7","3","non-coding","-","-",.,.,.,.,.,.,.,"1516",.,.,.,"MT"
   8281,"CCCCCTCTAG","del","MT-NC7","1","non-coding","-","-",.,.,.,.,.,.,.,"3",.,.,.,"MT"
   Update the reference count in the 8281 record, then delete the 8280 record
   perl -pi -e 's/8281,"CCCCCTCTA","del","MT-NC7","3"/8281,"CCCCCTCTA","del","MT-NC7","84"/' mitylib/annot/mitomap_panel_annotations.csv
   fgrep -v '8280,"CCCCCTCTA","del","MT-NC7","84",' mitylib/annot/mitomap_panel_annotations.csv > a; mv a mitylib/annot/mitomap_panel_annotations.csv


* re-run the above R code + vt normalize
These are the remaining errors... looks like i have not fixed some of the longer ones.

reference bases not consistent: MT:513-515  CCA(REF) vs CAC(FASTA)
reference bases not consistent: MT:514-516  AAC(REF) vs ACA(FASTA)
reference bases not consistent: MT:2351-2353  TTA(REF) vs TAA(FASTA)
reference bases not consistent: MT:3105-3106  CC(REF) vs CN(FASTA)
reference bases not consistent: MT:5133-5135  AAA(REF) vs ACT(FASTA)
reference bases not consistent: MT:8277-8286  TCTACCCCCT(REF) vs CTACCCCCTC(FASTA)
reference bases not consistent: MT:8287-8296  TCCCCCTCTA(REF) vs TAGAGCCCAC(FASTA)
reference bases not consistent: MT:9479-9494  TTTTTTCTTCGCAGGA(REF) vs TTTTTCTTCGCAGGAT(FASTA)
reference bases not consistent: MT:11620-11622  TTA(REF) vs TAC(FASTA)
reference bases not consistent: MT:15257-15258  GA(REF) vs AC(FASTA)
reference bases not consistent: MT:15496-15519  GCGACCCAGACAATTATACCCTAG(REF) vs GGCGACCCAGACAATTATACCCTA(FASTA)

* The FASTA records do not look right for some records that i checked. do this enmasse:
faidx mitylib/reference/hs37d5.MT.fa MT:513-515     | grep -v MT = GCA              -vs- CCA(REF) vs CAC(FASTA)
faidx mitylib/reference/hs37d5.MT.fa MT:514-516     | grep -v MT = CAC              -vs- AAC(REF) vs ACA(FASTA)
faidx mitylib/reference/hs37d5.MT.fa MT:2351-2353   | grep -v MT = TTA              -vs- TTA(REF) vs TAA(FASTA); "2352 TA-del"
faidx mitylib/reference/hs37d5.MT.fa MT:3105-3106   | grep -v MT = AC               -vs- CC(REF) vs CN(FASTA)
faidx mitylib/reference/hs37d5.MT.fa MT:5133-5135   | grep -v MT = AAC              -vs- AAA(REF) vs ACT(FASTA)
faidx mitylib/reference/hs37d5.MT.fa MT:8277-8286   | grep -v MT = TCTACCCCCT       -vs- TCTACCCCCT(REF) vs CTACCCCCTC(FASTA)
faidx mitylib/reference/hs37d5.MT.fa MT:8287-8296   | grep -v MT = CTAGAGCCCA       -vs- TCCCCCTCTA(REF) vs TAGAGCCCAC(FASTA)
faidx mitylib/reference/hs37d5.MT.fa MT:9479-9494   | grep -v MT = TTTTTTCTTCGCAGGA -vs- TTTTTTCTTCGCAGGA(REF) vs TTTTTCTTCGCAGGAT(FASTA)
faidx mitylib/reference/hs37d5.MT.fa MT:11620-11622 | grep -v MT = ATA              -vs- TTA(REF) vs TAC(FASTA)
faidx mitylib/reference/hs37d5.MT.fa MT:15257-15258 | grep -v MT = GA               -vs- GA(REF) vs AC(FASTA)
faidx mitylib/reference/hs37d5.MT.fa MT:15496-15519 | grep -v MT = GCA              -vs- GCGACCCAGACAATTATACCCTAG(REF) vs GGCGACCCAGACAATTATACCCTA(FASTA)

* This is proving very timeconsuming to fix records which are essentially wrong in Mitomap.
* we've fixed ~400 del records, so its time to move on.

# normalize & sort the variants using vt & linking back to R data.frame via the ID field
vcf <- with(orig,
    data.frame(
        CHROM="MT", 
        POS=POS,
        ID=1:nrow(orig),
        REF=REF,
        ALT=ALT,
        QUAL=".",
        FILTER=".",
        INFO=".")
    )
write.delim(vcf, "mitomap_panel_annotations_fixed.vcf")
# insert:
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=20191107
##reference=/Users/mcowley/src/mity/mitylib/reference/hs37d5.MT.fa							
##contig=<ID=MT,length=16569>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO


    $ vt normalize -n -r mitylib/reference/hs37d5.MT.fa -o mitomap_panel_annotations_fixed.left.vcf mitomap_panel_annotations_fixed.vcf 2>&1  | grep -c skipped
    ...
    no. right trimmed and left aligned    : 1
    no. left aligned                      : 382

    $ vt sort mitomap_panel_annotations_fixed.left.vcf -o mitomap_panel_annotations_fixed.left.sorted.vcf

    $ vt uniq mitomap_panel_annotations_fixed.left.sorted.vcf | wc -l
    uniq v0.57

    options:     input VCF file        mitomap_panel_annotations_fixed.left.sorted.vcf
             [o] output VCF file       -

    stats: Total number of observed variants   12150
           Total number of unique variants     12035

    Time elapsed: 0.02s

       12042
This wants to drop 115 variants. which record will be retained though by vt? It may drop a record with lots of disease support in mitomap.

vcf <- read.delim("mitomap_panel_annotations_fixed.left.sorted.vcf", skip=6)
colnames(vcf)[1] <- "CHROM"
orig$ID <- 1:nrow(orig)

res <- merge(
    vcf[,c("CHROM", "POS", "ID", "REF", "ALT", "INFO")], 
    orig[,c("ID", "locus_mitomap", "num_references_mitomap", 
            "variant_amino_acid_change_mitomap", "codon_position_mitomap", 
            "codon_number_mitomap", "disease_mitomap", "num_disease_references_mitomap", 
            "RNA_mitomap", "homoplasmy_mitomap", "heteroplasmy_mitomap", 
            "status_mitomap", "disease_amino_acid_change_mitomap", "GenBank_frequency_mitomap", 
            "baylor_panel", "common_22_panel", "common_58_panel",
            "commercial_panels")
        ]
    , by = "ID"
)

# tag dups; the x | y should tag the first and last record as DUP=TRUE which will work well for pairs. 
res$DUP <- duplicated(res[,c("POS", "REF", "ALT")], fromLast=TRUE) | duplicated(res[,c("POS", "REF", "ALT")], fromLast=FALSE)

#It's a tricky algo to pick the best record for a dup. Take the first variant for example:
     ID CHROM POS     REF   ALT                         INFO         locus_mitomap num_references_mitomap variant_amino_acid_change_mitomap
138 138    MT  65      TG     T                            .                MT-HV2                      1                                 .
150 150    MT  65      TG     T       OLD_VARIANT=MT:70:GG/G                MT-HV2                      5                                 .
    codon_position_mitomap codon_number_mitomap disease_mitomap num_disease_references_mitomap RNA_mitomap homoplasmy_mitomap heteroplasmy_mitomap
138                      .                    .               .                              .           .                  .                    .
150                      .                    .               .                              .           .                  .                    .
    status_mitomap disease_amino_acid_change_mitomap GenBank_frequency_mitomap baylor_panel common_22_panel common_58_panel commercial_panels  DUP
138              .                                 .                        20            .               .               .                  TRUE
150              .                                 .                         0            .               .               .                  TRUE

# The first row has higher GenBank_frequency_mitomap and lower num_references_mitomap...

# How many have been associated with a disease?
table(res[res$DUP,"status_mitomap"])
               .         Reported   Reported;.;.;. Reported;.;.;.;.          Unclear 
             209                2                2                1                1 
res[res$DUP & strlen(res$status_mitomap) > 3,]
         ID CHROM   POS REF              ALT                                    INFO                      locus_mitomap num_references_mitomap
568     568    MT   302   A              ACC                                       . MT-CR.MT-CSB2.MT-HV2.MT-OHR.MT-TFY              .;2.2.2.2
600     600    MT   302   A               AC                 OLD_VARIANT=MT:309:C/CC        MT-CR.MT-CSB2.MT-HV2.MT-OHR             .;71.71.71
629     629    MT   310   T               TC                                       .        MT-CR.MT-CSB2.MT-HV2.MT-OHR             .;20.20.20
1484   1484    MT   960  CT                C                                       .                            MT-RNR1                      .
10957 10957    MT 16015   T TATTCTCTGTTCTTTC OLD_VARIANT=MT:16032:T/TTCTCTGTTCTTTCAT                              MT-TP                      .
10959 10959    MT 16015   T TATTCTCTGTTCTTTC OLD_VARIANT=MT:16033:G/TCTCTGTTCTTTCATG                              MT-TP                      .
      variant_amino_acid_change_mitomap codon_position_mitomap codon_number_mitomap                                        disease_mitomap
568                                   .                      .                    .               Higher in melanoma patient group;.;.;.;.
600                                   .                      .                    .                             AD-weakly associated;.;.;.
629                                   .                      .                    .                                Melanoma patients;.;.;.
1484                                  .                      .                    .                                   DEAF / AD-associated
10957                                 .                      .                    .                     Dilated cardiomyopathy (15 bp dup)
10959                                 .                      .                    . Dilated cardiomyopathy (15 bp dup), alternate notation
      num_disease_references_mitomap RNA_mitomap homoplasmy_mitomap heteroplasmy_mitomap   status_mitomap disease_amino_acid_change_mitomap
568                        1;.;.;.;.           .                  .                    . Reported;.;.;.;.                 noncoding;.;.;.;.
600                          1;.;.;.           .                  .                    .   Reported;.;.;.                   noncoding;.;.;.
629                          1;.;.;.           .                  .                    .   Reported;.;.;.                   noncoding;.;.;.
1484                              18     12S rR.                  +                    +          Unclear                                 .
10957                              1     tR. Pro                  -                    +         Reported                                 .
10959                              1     tR. Pro                  -                    +         Reported                                 .
      GenBank_frequency_mitomap baylor_panel common_22_panel common_58_panel commercial_panels  DUP
568                          66            .               .               .                  TRUE
600                         284            .               .               .                  TRUE
629                           0            .               .               .                  TRUE
1484                          0            .               .               .                  TRUE
10957                         1            .               .               .                  TRUE
10959                         0            .               .               .                  TRUE
# wait, these ALT's are different; bow can the be tagged as DUPs?
# use vt unique and be done with it!

$ vt uniq mitomap_panel_annotations_fixed.left.sorted.vcf -o mitomap_panel_annotations_fixed.left.sorted.uniq.vcf
vcfu <- read.delim("mitomap_panel_annotations_fixed.left.sorted.uniq.vcf", skip=6)
colnames(vcfu)[1] <- "CHR"

res2 <- merge(
    vcfu[,c("CHR", "POS", "ID", "REF", "ALT", "INFO")], 
    orig[,c("ID", "locus_mitomap", "num_references_mitomap", 
            "variant_amino_acid_change_mitomap", "codon_position_mitomap", 
            "codon_number_mitomap", "disease_mitomap", "num_disease_references_mitomap", 
            "RNA_mitomap", "homoplasmy_mitomap", "heteroplasmy_mitomap", 
            "status_mitomap", "disease_amino_acid_change_mitomap", "GenBank_frequency_mitomap", 
            "baylor_panel", "common_22_panel", "common_58_panel",
            "commercial_panels")
        ]
    , by = "ID", all.x=T, all.y=F, no.dups=T
)

# reorder columns incase mity report deponds on the original order
res2 <- res2[, setdiff(colnames(orig), "ID")]
write.csv(res2, "dev/mitomap_panel_annotations_fixed.csv", row.names=FALSE)

# overwrite the original!
cp dev/mitomap_panel_annotations_fixed.csv  mitylib/annot/mitomap_panel_annotations.csv