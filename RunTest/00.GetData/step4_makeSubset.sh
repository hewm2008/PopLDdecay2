#!/bin/sh
# step4: rebuild the 50k slice that ships with this package (normally NOT needed)
#
#   sh step4_makeSubset.sh chr22/chr22.vcf.gz 50000
#
# Why the slice is filtered first: chr22 has 1,103,547 records, but only 217,103
# of them are SNPs that survive PopLDdecay2's own filters (44,030 indels and
# 7,835 multi-allelic sites are skipped, then -MAF 0.005 removes about 80% of the
# rest). Sampling the raw records would give a slice whose real site count is a
# fifth of its file name. Sampling the passing set instead gives a pure
# bi-allelic SNP slice with no skipped sites, so the site count is exact.
#
#   1. let PopLDdecay2 pick the sites that pass -MAF 0.005 / -Miss 0.25
#   2. take every step-th of them
#   3. pull those positions back out of the full VCF, with genotypes
#
# Result: 50,000 bi-allelic SNPs x 2,504 samples, spanning the whole chromosome.

VCF=${1:-chr22/chr22.vcf.gz}
M=${2:-50000}

W=$(mktemp -d)
trap 'rm -rf $W' EXIT

echo Start Time :
date

echo "1/4 sites passing -MAF 0.005 / -Miss 0.25"
../../bin/PopLDdecay2 -InVCF $VCF -OutStat $W/ref -MaxDist 300 -T 16 -OutFilterSNP >/dev/null 2>&1
zcat $W/ref.vcf.filter.gz | grep -v '^#' | cut -f2 > $W/pass.pos
echo "    $(wc -l < $W/pass.pos)"

echo "2/4 pick every step-th site"
awk -v m=$M -v n=$(wc -l < $W/pass.pos) 'NR % int(n/m) == 1' $W/pass.pos > $W/cand.pos
echo "    candidates $(wc -l < $W/cand.pos)"

echo "3/4 take them back out of the VCF, bi-allelic SNPs only"
zcat $VCF | awk -v m=$M '
	NR==FNR { want[$1]=1; next }
	/^#/ { print; next }
	($2 in want) && length($4)==1 && length($5)==1 &&
	$4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/ && $4 != $5 {
		print; if (++n >= m) exit
	}' $W/cand.pos - | bgzip -l 9 -c > base_chr22_50k.vcf.gz
tabix -f -p vcf base_chr22_50k.vcf.gz

echo "4/4 done"
ls -l base_chr22_50k.vcf.gz
echo "sites : $(zcat base_chr22_50k.vcf.gz | grep -vc '^#')"
echo "md5   : $(md5sum base_chr22_50k.vcf.gz | cut -d' ' -f1)"

echo End Time :
date
