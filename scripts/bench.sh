#!/bin/bash
# Benchmark harness with two explicit modes (public 1000 Genomes chr22 data):
#
#   bash scripts/bench.sh accuracy   (default)  SMALL data -> correctness gate
#     Old PopLDdecay vs PopLDdecay2 on a head extract of
#     Test/1000genomes_chr22/*.vcf.gz (2,504 samples, first N sites); every -T
#     value must be byte-identical (gzip -cd diff == 0). Small data verifies
#     accuracy quickly (single-chromosome, big-N=2504 path).
#
#   bash scripts/bench.sh perf                  BIG data -> performance
#     PopLDdecay2 only on the full Test/1000genomes_chr22/*.vcf.gz (chr22, 2,504
#     samples, 1,103,547 SNP); the old tool is very slow here (55 min), so the
#     baseline is the new tool itself at -T 1. Records wall/RSS per -T, checks
#     -T>1 vs -T 1 byte-consistency, reports throughput + scaling efficiency.
#     Single-chromosome input exercises the fastpath position-shard path (0.6.5),
#     so THREADS="1 8 16 32" covers sharded runs.
#     (NOT part of the per-build regression.)
#
# Overrides: VCF=..., THREADS="1 8 16 32"
set -uo pipefail
cd "$(dirname "$0")/.."

OLD=PopLDdecay-3.45/bin/PopLDdecay
NEW=bin/PopLDdecay2
WORK=Test//ldbench
mkdir -p "$WORK"
BGZIP=/home/heweiming/01.Software/samtools-1.23/bin/bgzip
TABIX=/home/heweiming/01.Software/samtools-1.23/bin/tabix
SRC=${SRC:-Test/1000genomes_chr22/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz}

MODE=${1:-accuracy}
THREADS=${THREADS:-"1 8 16 32"}

case "$MODE" in
accuracy)
  [ -f "$SRC" ] || { echo "accuracy mode needs $SRC"; exit 1; }
  HEAD="$WORK/chr22.head.vcf.gz"
  if [ ! -f "$HEAD" ]; then
    gzip -cd "$SRC" | awk 'BEGIN{k=0} /^#/{print; next} {k++; if(k<=2000) print; else exit}' > "$WORK/chr22.head.vcf"
    "$BGZIP" -f -c "$WORK/chr22.head.vcf" > "$HEAD" && rm -f "$WORK/chr22.head.vcf"
    "$TABIX" -p vcf "$HEAD"
  fi
  VCF=$HEAD
  echo "== accuracy gate ($VCF: 2000 sites x 2504 samples, chr22 head) =="
  ok=1
  for T in $THREADS; do
    "$OLD" -InVCF "$VCF" -OutStat "$WORK/acc.old.t$T.stat" 2>/dev/null
    "$NEW" -InVCF "$VCF" -OutStat "$WORK/acc.new.t$T.stat" -T "$T" 2>/dev/null
    d=$(diff <(gzip -cd "$WORK/acc.old.t$T.stat.gz") <(gzip -cd "$WORK/acc.new.t$T.stat.gz") | wc -l)
    echo "  -T $T: diff_vs_old=$d"
    [ "$d" = "0" ] || ok=0
  done
  # OutType 3 also writes .LD.gz; check both files vs old at -T 1 (fallback)
  # and -T 4 (indexed fastpath), exercising the sharded .LD.gz part merge.
  for T in 1 4; do
    "$OLD" -InVCF "$VCF" -OutStat "$WORK/acc.old.ot3.t$T.stat" -OutType 3 2>/dev/null
    "$NEW" -InVCF "$VCF" -OutStat "$WORK/acc.new.ot3.t$T.stat" -OutType 3 -T "$T" 2>/dev/null
    d=$(diff <(gzip -cd "$WORK/acc.old.ot3.t$T.stat.gz") <(gzip -cd "$WORK/acc.new.ot3.t$T.stat.gz") | wc -l)
    dl=$(diff <(gzip -cd "$WORK/acc.old.ot3.t$T.LD.gz") <(gzip -cd "$WORK/acc.new.ot3.t$T.LD.gz") | wc -l)
    echo "  -T $T (OutType 3): diff_vs_old=$d .LD.gz=$dl"
    [ "$d" = "0" ] && [ "$dl" = "0" ] || ok=0
  done
  [ "$ok" = "1" ] && echo "accuracy gate PASS" || echo "accuracy gate FAIL"
  exit $((1 - ok))
  ;;
perf)
  VCF=${VCF:-$SRC}
  [ -f "$VCF" ] || { echo "perf mode needs $VCF"; exit 1; }
  wall_secs() { echo "$1" | awk -F: '{ if (NF==3) print $1*3600+$2*60+$3; else print $1*60+$2 }'; }
  echo "== perf (self-baseline: PopLDdecay2 -T 1) =="
  echo "dataset: $VCF"
  # reference -T 1 run (baseline for byte-consistency)
  /usr/bin/time -v "$NEW" -InVCF "$VCF" -OutStat "$WORK/perf.t1.stat" -T 1 2> "$WORK/perf.t1.time" >/dev/null
  w1=$(wall_secs "$(grep 'Elapsed' "$WORK/perf.t1.time" | awk '{print $8}')")
  r1=$(grep 'Maximum resident' "$WORK/perf.t1.time" | awk '{print $6}')
  echo "  -T  1: wall=${w1}s rss=${r1}KB (baseline)"
  for T in $THREADS; do
    [ "$T" = "1" ] && continue
    /usr/bin/time -v "$NEW" -InVCF "$VCF" -OutStat "$WORK/perf.t$T.stat" -T "$T" 2> "$WORK/perf.t$T.time" >/dev/null
    w=$(wall_secs "$(grep 'Elapsed' "$WORK/perf.t$T.time" | awk '{print $8}')")
    r=$(grep 'Maximum resident' "$WORK/perf.t$T.time" | awk '{print $6}')
    v=$(diff <(gzip -cd "$WORK/perf.t1.stat.gz") <(gzip -cd "$WORK/perf.t$T.stat.gz") | wc -l)
    sp=$(awk -v a="$w1" -v b="$w" 'BEGIN{ if (b>0) printf "%.2fx", a/b; else print "n/a" }')
    eff=$(awk -v a="$w1" -v b="$w" -v t="$T" 'BEGIN{ if (b>0) printf "%.2f", a/(t*b); else print "n/a" }')
    echo "  -T $T: wall=${w}s rss=${r}KB speedup=${sp} eff=${eff} vs_T1=$v"
  done
  echo "  (throughput: see doc/BENCHMARK.md; run this on-demand, not per build)"
  exit 0
  ;;
*)
  echo "usage: bench.sh [accuracy|perf]"; exit 1
  ;;
esac