#!/bin/sh
# step2e: quickLD on every grid point (single-threaded)
#
# Its -threads option is non-functional in this release: it sets up the tasks
# and computes no pairs at all. quickLD also has no distance window on the
# correlation -- it emits every pair of the parsed variants -- so memory and
# output grow with the square of the site count; a cell whose projected memory
# does not fit is not started. The parse step (VCF -> .mdf) runs once per cell
# and is not timed; the timed step is qLD-compute.
#
# qLD-parse-VCF exits 1 even when it succeeds, so success is the header file
# it writes, not its status.

QUICKLD_DIR=${QUICKLD_DIR:-../../../rebuttal/tools/quickLD}
QD=$(cd "$QUICKLD_DIR" && pwd) || { echo "no quickLD in $QUICKLD_DIR -- set QUICKLD_DIR"; exit 1; }
SITES=${SITES:-"100000 200000 300000 500000"}
SAMPLES=${SAMPLES:-"250 500 1000 2504"}

mkdir -p out/quickld

echo Start Time :
date

for N in $SAMPLES; do
for M in $SITES; do
	V=$(pwd)/input/n${N}_m${M}.vcf.gz
	D=$(pwd)/out/quickld/n${N}_m${M}
	[ -f "$V" ] || { echo "no $V -- run step1_makeGrid.sh first"; exit 1; }
	need=$(awk -v n=$N -v m=$M 'BEGIN{printf "%.0f", 1584*n*m/1048576}')
	avail=$(awk '/MemAvailable/{printf "%.0f", $2/1024}' /proc/meminfo)
	[ $((need * 3)) -gt "$avail" ] && { echo "SKIPPED n${N}_m${M}: quickLD needs ~${need} MB, ${avail} MB free"; continue; }
	[ -s "$V.tbi" ] || tabix -f -p vcf "$V"
	mkdir -p $D
	if [ ! -d $D/qld.mdf ]; then
		echo "parse    n${N}_m${M}"
		( cd $D && rm -rf parsed parsed_1 \
			&& $QD/bin/qLD-parse-VCF -input $V -output parsed -chrom 22 -size 8192 > parse.log 2>&1
		  H=$(find . -name '*_header.vcf.gz' | head -1)
		  [ -n "$H" ] && $QD/bin/qLD-parse-2MDF -input "$(dirname $H)" -output qld.mdf >> parse.log 2>&1 )
		[ -d $D/qld.mdf ] || { echo "parse failed for n${N}_m${M} (see $D/parse.log)"; continue; }
	fi
	echo "quickLD  n${N}_m${M}"
	( cd $D && /bin/time	-v	-p	$QD/bin/qLD-compute	-input	qld.mdf	-output	qld.out	-ploidy	phased_diploid	-r2limit	0	-mdf	-blis	2>	qld.time )
done
done

echo End Time :
date
