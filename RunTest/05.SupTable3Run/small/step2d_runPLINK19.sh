#!/bin/sh
# step2d: PLINK 1.9 on every grid point, --r2 square
#
# Square mode holds a kept x kept matrix: 20,697 kept variants measured
# 3,339 MB, and the footprint grows with the square of the kept count. A cell
# that does not fit in three times the free memory is not started -- this is a
# shared node. The bed is the one step2c built; the conversion is not timed.

PLINK2=${PLINK2:-$(command -v plink2 || echo /home/heweiming/01.Software/plink/plink2)}
PLINK19=${PLINK19:-$(command -v plink || echo /home/heweiming/01.Software/plink-1.90b6.21/bin/plink)}
SITES=${SITES:-"2000 5000 10000 50000 100000"}
SAMPLES=${SAMPLES:-"250 500 1000 2504"}
REPEATS=${REPEATS:-5}
THREADS=${THREADS:-8}

mkdir -p out/bed out/plink19

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
	k=$(wc -l < $BED.bim)
	R=1
	while [ "$R" -le "$REPEATS" ]; do
		need=$(awk -v k=$k 'BEGIN{printf "%.0f", 3339*(k/20697)^2}')
		avail=$(awk '/MemAvailable/{printf "%.0f", $2/1024}' /proc/meminfo)
		if [ $((need * 3)) -gt "$avail" ]; then
			echo "SKIPPED n${N}_m${M} rep $R: PLINK 1.9 needs ~${need} MB, ${avail} MB free"
			printf 'SKIPPED_MEMORY>need_%sMB_avail_%sMB\n' "$need" "$avail" > out/plink19/n${N}_m${M}.r$R.status
			R=$((R + 1))
			continue
		fi
		echo "PLINK 1.9  n${N}_m${M}  ($k kept variants)  rep $R/$REPEATS"
		/bin/time	-v	-p	$PLINK19	--bfile	$BED	--r2	square	gz	--threads	$THREADS	--out	out/plink19/n${N}_m${M}	>	out/plink19/n${N}_m${M}.r$R.log	2>	out/plink19/n${N}_m${M}.r$R.time
		rm -f out/plink19/n${N}_m${M}.ld.gz
		R=$((R + 1))
	done
done
done

echo End Time :
date
