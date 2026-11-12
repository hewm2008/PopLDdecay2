#!/bin/sh
# make the 16 "large" grid datasets
#
#   samples : 250 / 500 / 1000 / 2504
#   sites   : 100,000 / 200,000 / 300,000 / 500,000
#
# Every cell is a systematic sample of the biallelic SNPs of chr22: record i
# of the T biallelic SNPs is kept when (i*M) mod T < M, which is true exactly
# M times, so the M sites are spread over the whole chromosome instead of
# being taken from its start, and the spacing between neighbouring sites is
# preserved. The first N sample columns are then cut out.
#
# Only biallelic SNPs are used: PLINK 2 refuses a .bim with multiallelic
# variants, and PopLDdecay2, PopLDdecay 3.45 and LDkit skip them, so a grid
# with multiallelic sites would not give the tools the same input. These 16
# cells are the grid of the published Table S3.
#
# Needs chr22 (../../00.GetData/step1_download.sh) and its rebuilt index.
# The 16 files total roughly 1 GB.

CHR22=${CHR22:-../../00.GetData/chr22/chr22.vcf.gz}
SITES=${SITES:-"100000 200000 300000 500000"}
SAMPLES=${SAMPLES:-"250 500 1000 2504"}
BGZIP=${BGZIP:-$(command -v bgzip || echo /home/heweiming/01.Software/samtools-1.23/bin/bgzip)}
TABIX=${TABIX:-$(command -v tabix || echo /home/heweiming/01.Software/samtools-1.23/bin/tabix)}

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
		| $BGZIP -c > input/n${N}_m${M}.vcf.gz
	$TABIX -@ 8 -f -p vcf input/n${N}_m${M}.vcf.gz
	echo "n${N}_m${M}  sites=$(zcat input/n${N}_m${M}.vcf.gz | grep -vc '^#')  samples=$(zcat input/n${N}_m${M}.vcf.gz | grep -m1 '^#CHROM' | awk '{print NF-9}')"
done
done

echo End Time :
date
echo "next: sh step2_runGrid.sh"
