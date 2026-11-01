#!/bin/sh
# step1: download 1000 Genomes Phase 3 (20130502) chr22 and the sample panel
#
# Two sources, same content:
#   EBI  (v5b, used for the paper)
#     https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
#   UCSC (v5a only -- a DIFFERENT site set, so the numbers will not match)
#     https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/
#
# Already have a local copy? Pass its directory:
#   sh step1_download.sh /path/to/1000genomes_20130502

echo Start Time :
date

mkdir -p chr22

if [ $# -ge 1 ]; then
	cp $1/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz chr22/chr22.vcf.gz
	cp $1/integrated_call_samples_v3.20130502.ALL.panel chr22/panel.txt
else
	wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
		-O chr22/chr22.vcf.gz
	wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel \
		-O chr22/panel.txt
fi

echo "sites   : $(zcat chr22/chr22.vcf.gz | grep -vc '^#')"
echo "samples : $(zcat chr22/chr22.vcf.gz | grep -m1 '^#CHROM' | awk '{print NF-9}')"

echo End Time :
date
echo "next: sh step3_makeTBI.sh"
