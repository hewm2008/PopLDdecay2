#!/bin/sh
# step2b: PopLDdecay 3.45 on every grid point (single-threaded by design)
#
# The original tool has no -T and no indexed fast path. Same input, 300 kb
# window and MAF 0.005 as the other six tools.

PLD345=../../../PopLDdecay-3.45/bin/PopLDdecay
SITES=${SITES:-"2000 5000 10000 50000 100000"}
SAMPLES=${SAMPLES:-"250 500 1000 2504"}
REPEATS=${REPEATS:-5}

mkdir -p out/pld345

echo Start Time :
date

for N in $SAMPLES; do
for M in $SITES; do
	V=input/n${N}_m${M}.vcf.gz
	[ -f "$V" ] || { echo "no $V -- run step1_makeGrid.sh first"; exit 1; }
	R=1
	while [ "$R" -le "$REPEATS" ]; do
		echo "PopLDdecay 3.45  n${N}_m${M}  rep $R/$REPEATS"
		/bin/time	-v	-p	$PLD345	-InVCF	$V	-OutStat	out/pld345/n${N}_m${M}.old	-MaxDist	300	-MAF	0.005	>	out/pld345/n${N}_m${M}.r$R.log	2>	out/pld345/n${N}_m${M}.r$R.time
		R=$((R + 1))
	done
done
done

echo End Time :
date
