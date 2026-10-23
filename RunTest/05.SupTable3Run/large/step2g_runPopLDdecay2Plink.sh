#!/bin/sh
# step2g: PopLDdecay2 reading PLINK input (-InPlink), -T 1 and -T 8
#
# Runs on the same bed as PLINK 2 and PLINK 1.9 (step2c), so the engines read
# one identical binary dataset and neither pays a format cost -- the matched
# comparison in which the VCF parsing of step2a is out of the picture.
# -InPlink has no indexed fast path, so -T only parallelises the decoding; and
# a .bed is unphased, so this step computes the unphased r^2.

PLD2=../../../bin/PopLDdecay2
PLINK2=${PLINK2:-$(command -v plink2 || echo /home/liaojing/01.Software/plink/plink2)}
SITES=${SITES:-"100000 200000 300000 500000"}
SAMPLES=${SAMPLES:-"250 500 1000 2504"}

mkdir -p out/bed out/pld2plink

echo Start Time :
date

for N in $SAMPLES; do
for M in $SITES; do
	V=input/n${N}_m${M}.vcf.gz
	BED=out/bed/n${N}_m${M}
	[ -f "$V" ] || { echo "no $V -- run step1_makeGrid.sh first"; exit 1; }
	if [ ! -f "$BED.bed" ]; then
		echo "convert  n${N}_m${M}"
		/bin/time	-v	-p	$PLINK2	--vcf-half-call	missing	--make-bed	--vcf	$V	--out	$BED	--maf	0.005	--geno	0.25	2>	$BED.convert.time
	fi
	for T in 1 8; do
		echo "PopLDdecay2 -InPlink  n${N}_m${M}  -T $T"
		/bin/time	-v	-p	$PLD2	-InPlink	$BED.bed	-OutStat	out/pld2plink/n${N}_m${M}.t$T	-MaxDist	300	-MAF	0.005	-T	$T	2>	out/pld2plink/n${N}_m${M}.t$T.time
	done
done
done

echo End Time :
date
