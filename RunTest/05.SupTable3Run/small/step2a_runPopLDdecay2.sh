#!/bin/sh
# step2a: PopLDdecay2 on every grid point, -T 1 and -T 8
#
# -T 1 is the matched single-threaded column of Table S3 and -T 8 the
# multi-threaded one. Same input/n<N>_m<M>.vcf.gz, 300 kb window and MAF 0.005
# as the other six tools. The .tbi index is what makes -T 8 fast; without it
# PopLDdecay2 falls back to a sequential read.

PLD2=../../../bin/PopLDdecay2
SITES=${SITES:-"2000 5000 10000 50000 100000"}
SAMPLES=${SAMPLES:-"250 500 1000 2504"}

mkdir -p out/pld2

echo Start Time :
date

for N in $SAMPLES; do
for M in $SITES; do
	V=input/n${N}_m${M}.vcf.gz
	[ -f "$V" ] || { echo "no $V -- run step1_makeGrid.sh first"; exit 1; }
	[ -s "$V.tbi" ] || tabix -f -p vcf "$V"
	for T in 1 8; do
		echo "PopLDdecay2  n${N}_m${M}  -T $T"
		/bin/time	-v	-p	$PLD2	-InVCF	$V	-OutStat	out/pld2/n${N}_m${M}.t$T	-MaxDist	300	-MAF	0.005	-T	$T	2>	out/pld2/n${N}_m${M}.t$T.time
	done
done
done

echo End Time :
date
