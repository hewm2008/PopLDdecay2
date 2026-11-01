#!/bin/sh
# step2e: quickLD on every grid point. One MDF task per cell, so -competing
# avoids the chunk=0 scheduler bug. No distance window: memory/output grow with
# sites^2. Parse is not timed. qLD-parse-VCF exits 1 on success; check header.

QUICKLD_DIR=${QUICKLD_DIR:-../../../rebuttal/tools/quickLD}
QD=$(cd "$QUICKLD_DIR" && pwd) || { echo "no quickLD in $QUICKLD_DIR -- set QUICKLD_DIR"; exit 1; }
SITES=${SITES:-"2000 5000 10000 50000 100000"}
SAMPLES=${SAMPLES:-"250 500 1000 2504"}
REPEATS=${REPEATS:-5}
THREADS=${THREADS:-8}
QLD_TIMEOUT=${QLD_TIMEOUT:-14400}
TABIX=${TABIX:-$(command -v tabix || echo /home/heweiming/01.Software/samtools-1.23/bin/tabix)}

write_status() {
	R=1
	while [ "$R" -le "$REPEATS" ]; do echo "$1" > "$D/qld.r$R.status"; R=$((R + 1)); done
}

skip_cell() {
	echo "$1"
	write_status "$2"
}

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
		skip_cell "SKIPPED n${N}_m${M}: quickLD not run above 200k sites" "SKIPPED_SIZE>${M}_sites"
		continue
		;;
	esac
	need=$(awk -v n=$N -v m=$M 'BEGIN{printf "%.0f", 1584*n*m/1048576}')
	avail=$(awk '/MemAvailable/{printf "%.0f", $2/1024}' /proc/meminfo)
	if [ $((need * 3)) -gt "$avail" ]; then
		skip_cell "SKIPPED n${N}_m${M}: quickLD needs ~${need} MB, ${avail} MB free" "SKIPPED_MEMORY>need_${need}MB_avail_${avail}MB"
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
			write_status PARSE_FAILED
			continue
		fi
	fi
	R=1
	while [ "$R" -le "$REPEATS" ]; do
		echo "quickLD  n${N}_m${M}  rep $R/$REPEATS"
		rm -f $D/qld.out $D/qld.r$R.status
		( cd $D && /bin/time	-v	-p	timeout	"$QLD_TIMEOUT"	$QD/bin/qLD-compute	-input	qld.mdf	-output	qld.out	-ploidy	phased_diploid	-r2limit	0	-mdf	-blis	-threads	$THREADS	-competing	> qld.r$R.log 2> qld.r$R.time )
		calc=$(awk '/lds calculated/{for(i=1;i<=NF;i++) if($i=="lds") s+=$(i-1)} END{printf "%.0f", s+0}' $D/qld.r$R.log 2>/dev/null)
		rc=$(sed -n 's/.*Exit status: *\([0-9][0-9]*\).*/\1/p' $D/qld.r$R.time | tail -1)
		status=
		case "$rc" in
		124) status="TIMEOUT>${QLD_TIMEOUT}s" ;;
		'') status='FAILED<no_exit_status>' ;;
		0) ;;
		*) if [ "$rc" -gt 128 ]; then status="CRASH_SIGNAL$((rc - 128))"; else status="FAILED<rc${rc}>"; fi ;;
		esac
		[ -z "$status" ] && { [ ! -s $D/qld.out ] || [ "${calc:-0}" -eq 0 ]; } && status="EMPTY_OUTPUT>lds_${calc:-0}"
		[ -n "$status" ] && echo "$status" > "$D/qld.r$R.status"
		R=$((R + 1))
	done
done
done

echo End Time :
date
