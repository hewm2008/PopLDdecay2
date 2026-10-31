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
REPEATS=${REPEATS:-5}
THREADS=${THREADS:-8}
QLD_TIMEOUT=${QLD_TIMEOUT:-14400}
TABIX=${TABIX:-$(command -v tabix || echo /home/heweiming/01.Software/samtools-1.23/bin/tabix)}

mkdir -p out/quickld

echo Start Time :
date

for N in $SAMPLES; do
for M in $SITES; do
	V=$(pwd)/input/n${N}_m${M}.vcf.gz
	D=$(pwd)/out/quickld/n${N}_m${M}
	[ -f "$V" ] || { echo "no $V -- run step1_makeGrid.sh first"; exit 1; }
	mkdir -p $D
	case "$M" in
	300000|500000)
		echo "SKIPPED n${N}_m${M}: quickLD not run above 200k sites"
		R=1
		while [ "$R" -le "$REPEATS" ]; do
			echo "SKIPPED_SIZE>${M}_sites" > $D/qld.r$R.status
			R=$((R + 1))
		done
		continue
		;;
	esac
	need=$(awk -v n=$N -v m=$M 'BEGIN{printf "%.0f", 1584*n*m/1048576}')
	avail=$(awk '/MemAvailable/{printf "%.0f", $2/1024}' /proc/meminfo)
	if [ $((need * 3)) -gt "$avail" ]; then
		echo "SKIPPED n${N}_m${M}: quickLD needs ~${need} MB, ${avail} MB free"
		R=1
		while [ "$R" -le "$REPEATS" ]; do
			printf 'SKIPPED_MEMORY>need_%sMB_avail_%sMB\n' "$need" "$avail" > $D/qld.r$R.status
			R=$((R + 1))
		done
		continue
	fi
	[ -s "$V.tbi" ] || $TABIX -@ 8 -f -p vcf "$V"
	if [ ! -d $D/qld.mdf ]; then
		echo "parse    n${N}_m${M}"
		( cd $D && rm -rf parsed parsed_1 \
			&& timeout "$QLD_TIMEOUT" "$QD/bin/qLD-parse-VCF" -input $V -output parsed -chrom 22 -size 8192 > parse.log 2>&1
		  H=$(find . -name '*_header.vcf.gz' | head -1)
		  [ -n "$H" ] && timeout "$QLD_TIMEOUT" "$QD/bin/qLD-parse-2MDF" -input "$(dirname $H)" -output qld.mdf >> parse.log 2>&1 )
		if [ ! -d $D/qld.mdf ]; then
			echo "parse failed for n${N}_m${M} (see $D/parse.log)"
			R=1
			while [ "$R" -le "$REPEATS" ]; do
				echo PARSE_FAILED > $D/qld.r$R.status
				R=$((R + 1))
			done
			continue
		fi
	fi
	R=1
	while [ "$R" -le "$REPEATS" ]; do
		echo "quickLD  n${N}_m${M}  rep $R/$REPEATS"
		( cd $D && /bin/time	-v	-p	timeout	"$QLD_TIMEOUT"	$QD/bin/qLD-compute	-input	qld.mdf	-output	qld.out	-ploidy	phased_diploid	-r2limit	0	-mdf	-blis	-threads	$THREADS	> qld.r$R.log 2> qld.r$R.time )
		calc=$(awk '/lds calculated/{for(i=1;i<=NF;i++) if($i=="lds") s+=$(i-1)} END{printf "%.0f", s+0}' $D/qld.r$R.log 2>/dev/null)
		rc=$(sed -n 's/.*Exit status: *\([0-9][0-9]*\).*/\1/p' $D/qld.r$R.time | tail -1)
		if [ "$rc" != 124 ] && { [ ! -s $D/qld.out ] || [ "${calc:-0}" -eq 0 ]; }; then
			printf 'EMPTY_OUTPUT>lds_%s\n' "${calc:-0}" > $D/qld.r$R.status
		fi
		R=$((R + 1))
	done
done
done

echo End Time :
date
