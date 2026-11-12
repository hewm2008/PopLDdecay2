#!/bin/sh
# step2c: PLINK 2 on every grid point
#
# The VCF -> bed conversion is not part of the timed LD step (Table S3 reports
# the PLINK columns as a lower bound of the end-to-end cost); its record goes
# to out/bed/<cell>.convert.time. The bed is shared with step2d and step2g, so
# PLINK 1.9 and PopLDdecay2 -InPlink read the very same dataset.
#
# --r2-phased is the haplotype r^2, the statistic PopLDdecay2 computes.
# --ld-window-r2 0 is required: the default 0.2 would drop the low-r^2 pairs.

PLINK2=${PLINK2:-$(command -v plink2 || echo /home/heweiming/01.Software/plink/plink2)}
SITES=${SITES:-"2000 5000 10000 50000 100000"}
SAMPLES=${SAMPLES:-"250 500 1000 2504"}
REPEATS=${REPEATS:-5}
THREADS=${THREADS:-8}

mkdir -p out/bed out/plink2

echo Start Time :
date

for N in $SAMPLES; do
for M in $SITES; do
	V=input/n${N}_m${M}.vcf.gz
	BED=out/bed/n${N}_m${M}
	[ -f "$V" ] || { echo "no $V -- run step1_makeGrid.sh first"; exit 1; }
	if [ ! -f "$BED.bed" ]; then
		echo "convert  n${N}_m${M}"
		/bin/time	-v	-p	$PLINK2	--vcf-half-call	missing	--make-bed	--vcf	$V	--out	$BED	--maf	0.005	--geno	0.25	--threads	$THREADS	2>	$BED.convert.time
	fi
	R=1
	while [ "$R" -le "$REPEATS" ]; do
		echo "PLINK 2  n${N}_m${M}  rep $R/$REPEATS"
		/bin/time	-v	-p	$PLINK2	--bfile	$BED	--r2-phased	--ld-window-kb	300	--ld-window-r2	0	--threads	$THREADS	--out	out/plink2/n${N}_m${M}	>	out/plink2/n${N}_m${M}.r$R.log	2>	out/plink2/n${N}_m${M}.r$R.time
		R=$((R + 1))
	done
done
done

echo End Time :
date
