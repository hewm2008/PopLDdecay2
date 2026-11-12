#!/bin/sh
# step2f: LDkit on every grid point
#
# Same input, 300 kb window, MAF 0.005 and r^2 as the other tools. LDkit
# deletes the tabix index of the file it reads, which PopLDdecay2 and quickLD
# need, so the index is rebuilt after each run.

JAVA=${JAVA:-java}
LDKIT_JAR=${LDKIT_JAR:-../../../rebuttal/tools/LDkit.jar}
SITES=${SITES:-"2000 5000 10000 50000 100000"}
SAMPLES=${SAMPLES:-"250 500 1000 2504"}
REPEATS=${REPEATS:-5}
THREADS=${THREADS:-8}
TABIX=${TABIX:-$(command -v tabix || echo /home/heweiming/01.Software/samtools-1.23/bin/tabix)}

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
	R=1
	while [ "$R" -le "$REPEATS" ]; do
		rm -rf $D/ldkit
		echo "LDkit  n${N}_m${M}  rep $R/$REPEATS"
		/bin/time	-v	-p	$JAVA	-Xmx16g	-jar	$LDKIT_JAR	--infile	$V	--out	$D/ldkit	--ws	300	--maf	0.005	--type	1	--threads	$THREADS	> $D/ldkit.r$R.log 2> $D/ldkit.r$R.time
		[ -s "$V.tbi" ] || $TABIX -@ 8 -f -p vcf "$V"
		R=$((R + 1))
	done
done
done

echo End Time :
date
