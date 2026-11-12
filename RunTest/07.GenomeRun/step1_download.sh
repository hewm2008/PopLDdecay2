#!/bin/sh
# step1: fetch the 22 autosomes of 1000 Genomes Phase 3, build v5b (about 14 GB)
#
#   sh step1_download.sh              # from EBI
#   sh step1_download.sh /path/to/dir # from a local copy instead
#
# v5b is the build the paper used. The UCSC mirror carries only v5a, a
# different site set, so it will not reproduce these numbers.

P=phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
BASE=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502

mkdir -p genome

echo Start Time :
date

for i in $(seq 1 22); do
	F=ALL.chr${i}.$P
	if [ $# -ge 1 ]; then
		cp $1/$F genome/
	else
		wget -c -P genome $BASE/$F
	fi
done

ls -1 genome/*.vcf.gz | wc -l
du -sh genome

echo End Time :
date
