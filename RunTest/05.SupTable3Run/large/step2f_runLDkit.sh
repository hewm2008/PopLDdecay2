#!/bin/sh
# step2f: LDkit on every grid point
#
# Same input, 300 kb window, MAF 0.005 and r^2 as the other tools. LDkit
# deletes the tabix index of the file it reads, which PopLDdecay2 and quickLD
# need, so the index is rebuilt after each run.

JAVA=${JAVA:-java}
LDKIT_JAR=${LDKIT_JAR:-../../../rebuttal/tools/LDkit.jar}
SITES=${SITES:-"100000 200000 300000 500000"}
SAMPLES=${SAMPLES:-"250 500 1000 2504"}

mkdir -p out/ldkit

echo Start Time :
date

for N in $SAMPLES; do
for M in $SITES; do
	V=input/n${N}_m${M}.vcf.gz
	D=out/ldkit/n${N}_m${M}
	[ -f "$V" ] || { echo "no $V -- run step1_makeGrid.sh first"; exit 1; }
	[ -f "$LDKIT_JAR" ] || { echo "no $LDKIT_JAR -- set LDKIT_JAR"; exit 1; }
	mkdir -p $D
	echo "LDkit  n${N}_m${M}"
	/bin/time	-v	-p	$JAVA	-Xmx16g	-jar	$LDKIT_JAR	--infile	$V	--out	$D/ldkit	--ws	300	--maf	0.005	--type	1	--threads	8	2>	$D/ldkit.time
	[ -s "$V.tbi" ] || tabix -f -p vcf "$V"
done
done

echo End Time :
date
