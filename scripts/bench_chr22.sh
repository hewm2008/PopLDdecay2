#!/bin/bash
# Paper benchmark 1: chr22 1000 Genomes Phase 3 (2504 samples, ~1.05M SNP).
# PopLDdecay2 vs PopLDdecay-3.45: wall/RSS/output-size per run plus a
# decompressed .stat.gz byte-diff against the 3.45 reference. 3.45 is
# single-threaded (no -T); PopLDdecay2 is run at -T 1/4/8/16 to show both the
# single-thread speedup and thread scaling. Results land in Test//paper_chr22.
# Designed to run overnight in the background (setsid nohup ... &).
# NOTE: the chr22 VCF copied from the 1000 Genomes website ships with a STALE
# .tbi index (data rewritten after the index was built), which breaks the indexed
# fast path and makes PopLDdecay2 output differ from 3.45. Rebuild the index
# before running:
#   Test/1000genomes_chr22/*.vcf.gz.tbi   (delete the old .tbi)
#   /data/01.SoftWare/samtools-1.23/bin/tabix -f Test/1000genomes_chr22/*.vcf.gz
set -uo pipefail
cd "$(dirname "$0")/.."

VCF=${VCF:-Test/1000genomes_chr22/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz}
WORK=Test//paper_chr22
mkdir -p "$WORK"
OLD=PopLDdecay-3.45/bin/PopLDdecay
NEW=bin/PopLDdecay2
MAXD=${MAXD:-300}

wall_secs() { echo "$1" | awk -F: '{ if (NF==3) print $1*3600+$2*60+$3; else print $1*60+$2 }'; }

bench() {  # name bin extra...
  local name="$1" bin="$2"; shift 2
  /usr/bin/time -v "$bin" -InVCF "$VCF" -OutStat "$WORK/$name.stat" -MaxDist "$MAXD" "$@" \
    > "$WORK/$name.out" 2> "$WORK/$name.time"
  local rc=$? w r
  w=$(wall_secs "$(grep 'Elapsed' "$WORK/$name.time" | awk '{print $8}')")
  r=$(grep 'Maximum resident' "$WORK/$name.time" | awk '{print $6}')
  echo "  $name: rc=$rc wall=${w}s rss=$((r/1024))MB"
}

echo "== chr22 1000G benchmark (MaxDist=${MAXD}kb) =="
echo "data: $(zcat "$VCF" | grep -vc '^#') sites, $(zcat "$VCF" | grep -m1 '^#CHROM' | cut -f10- | wc -w) samples"
date +%s > "$WORK/start.ts"

bench old.345   "$OLD"
bench new.t1    "$NEW" -T 1
bench new.t4    "$NEW" -T 4
bench new.t8    "$NEW" -T 8
bench new.t16   "$NEW" -T 16

date +%s > "$WORK/end.ts"

echo "== result consistency (decompressed .stat.gz diff vs 3.45) =="
for name in new.t1 new.t4 new.t8 new.t16; do
  d=$(diff <(gzip -cd "$WORK/old.345.stat.gz") <(gzip -cd "$WORK/$name.stat.gz") | wc -l)
  echo "  $name vs 3.45: diff=$d"
done
echo "== output sizes =="
ls -l "$WORK"/*.stat.gz | awk '{print "  "$NF" "$5}'
echo "== elapsed total: $(( $(cat "$WORK/end.ts") - $(cat "$WORK/start.ts") ))s =="
