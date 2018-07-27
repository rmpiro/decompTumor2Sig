#!/bin/bash

# convert the MPF file to a VCF file
# use bcftools to convert the data for all samples

for PID in `zcat inst/extdata/Nik-Zainal_PMID_22608084-MPF.txt.gz | cut -f 1 | sort | uniq`; do 

	# do this for each patient/sample ID individually

	# create temporary file, adding a dummy header to be read by bcftools
	TSVFILE="tmp.$PID.tsv"; 
	echo -e "ID\tCHROM\tPOS\tAA" > $TSVFILE; 

	# get all mutations for the patientID, filter for SNVs
	zcat inst/extdata/Nik-Zainal_PMID_22608084-MPF.txt.gz | grep "^$PID" | grep "[[:space:]][ACGT][[:space:]][ACGT]$" | cut -f 2-5 | sed -e "s/[[:space:]]\([ACGT]\)[[:space:]]\([ACGT]\)$/\t\1\2/; s/^/.\t/; s/chr//" >> $TSVFILE; 

	# convert this TSV file to VCF using bcftools
	bcftools convert -c ID,CHROM,POS,AA -s $PID -f ~/bioinfo/data/ngs-data/ref-genomes/human/hs37d5.fa --tsv2vcf $TSVFILE -Oz -o $TSVFILE.vcf.gz ; 

	# remove the temporary TSV file
	rm -f $TSVFILE; 

	# index the temporary VCF file (tabix)
	tabix -p vcf $TSVFILE.vcf.gz; 

done; 

# now merge all the temporary VCF files into one multi-sample VCF file
vcf-merge tmp.*.tsv.vcf.gz | bgzip -c > Nik-Zainal_PMID_22608084-VCF-convertedfromMPF.vcf.gz; 

# remove the temporary VCF file for single samples
rm -f tmp.*.tsv.vcf.gz*


