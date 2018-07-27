#!/bin/bash

# Paper: https://www.nature.com/articles/nature12477
# Somatic mutations available at:
# ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/

# get all the somatic mutation data (about 767M) ...
wget -np -nH --cut-dirs=3 -r ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/


# get IDs for all genomes:
OUTF="genomes-tumor-patientID-map.tsv"; rm -f $OUTF; for TUMOR in `ls mutational_catalogs/genomes/ | sed -e "s/ .*//"`; do echo $TUMOR; cat mutational_catalogs/genomes/$TUMOR*/*subs.txt | grep "^Mutation Type" | sed -e "s/\t/\n/g;" | grep -v "^Mutation" | sed -e "s/^/$TUMOR\t/" >> $OUTF; done

cat genomes-tumor-patientID-map.tsv | cut -f 2 > genomes-patientID.list


# get all SNVs, but from genomes only:

cat somatic_mutation_data/*/*_clean_somatic_mutations_for_signature_analysis.txt | grep -w subs | ./extractSpecColumns.pl -c 3,1,4,6,7 | sed -e "s/^/chr/" | ./extractSpecColumns.pl -c 2,1,3-5 | filterLines.pl -c 1 -f genomes-patientID.list | gzip > Somatic_SNVs_Alexandrov_genomes_only.tsv.gz

# exclude MT:

zcat Somatic_SNVs_Alexandrov_genomes_only.tsv.gz | grep -w chrMT -v | gzip > Somatic_SNVs_Alexandrov_genomes_only-noMT.tsv.gz

# which samples have min 100 SNVs?
zcat Somatic_SNVs_Alexandrov_genomes_only-noMT.tsv.gz | cut -f 1 | sort | uniq -c | sed -e "s/^ *//; s/ /\t/" | sort -n | filterLinesNumeric.pl -c 1 -T GE -V 100 | cut -f 2 > genomes-patientID-min100snvs.list

# get the corresponding SNVs:

zcat Somatic_SNVs_Alexandrov_genomes_only-noMT.tsv.gz | filterLines.pl -c 1 -f genomes-patientID-min100snvs.list | gzip > Somatic_SNVs_Alexandrov_genomes_only-noMT-min100snvs.tsv.gz

# remove intermediate files (if you want):

rm -f genomes-tumor-patientID-map.tsv \
   genomes-patientID.list \
   Somatic_SNVs_Alexandrov_genomes_only.tsv.gz \
   Somatic_SNVs_Alexandrov_genomes_only-noMT.tsv.gz \
   genomes-patientID-min100snvs.list
   


# postprocess with pmsignature and decompTumor2Sig to get and save signatures

R CMD BATCH --vanilla Alexandrov_PMID_23945592_435_tumors-pmsignature-15sig.R
