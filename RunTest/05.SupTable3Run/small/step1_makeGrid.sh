#!/bin/sh
# make the 20 "small" grid datasets
#
#   samples : 250 / 500 / 1000 / 2504
#   sites   : 2,000 / 5,000 / 10,000 / 50,000 / 100,000
#
# Every cell is a systematic sample of the biallelic SNPs of chr22: record i
# of the T biallelic SNPs is kept when (i*M) mod T < M, which is true exactly
# M times, so the M sites are spread over the whole chromosome instead of
# being taken from its start, and the spacing between neighbouring sites is
# preserved. The first N sample columns are then cut out.
#
# Only biallelic SNPs are used: PLINK 2 refuses a .bim with multiallelic
# variants, and PopLDdecay2, PopLDdecay 3.45 and LDkit skip them, so a grid
# with multiallelic sites would not give the tools the same input.
#
# Needs chr22 (../../00.GetData/step1_download.sh) and its rebuilt index.

CHR22=${CHR22:-../../00.GetData/chr22/chr22.vcf.gz}
SITES=${SITES:-"2000 5000 10000 50000 100000"}
SAMPLES=${SAMPLES:-"250 500 1000 2504"}

mkdir -p input

echo Start Time :
date

# T = the number of biallelic SNPs on chr22 (1,055,454 for 1000G v5b)
T=$(zcat $CHR22 | awk 'BEGIN{FS="\t"} /^#/{next} {if(length($4)==1 && length($5)==1 && $4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/) n++} END{printf "%d", n}')
echo "biallelic SNPs on chr22: $T"

for N in $SAMPLES; do
for M in $SITES; do
	zcat $CHR22 \
		| awk -v m=$M -v t=$T 'BEGIN{FS="\t";i=0} /^#/{print;next} {if(length($4)==1 && length($5)==1 && $4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/){i++; if((i*m)%t < m) print}}' \
		| cut -f1-9,10-$((9 + N)) \
		| bgzip -c > input/n${N}_m${M}.vcf.gz
	tabix -f -p vcf input/n${N}_m${M}.vcf.gz
	echo "n${N}_m${M}  sites=$(zcat input/n${N}_m${M}.vcf.gz | grep -vc '^#')  samples=$(zcat input/n${N}_m${M}.vcf.gz | grep -m1 '^#CHROM' | awk '{print NF-9}')"
done
done

echo End Time :
date
echo "next: sh step2_runGrid.sh"
