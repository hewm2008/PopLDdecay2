#!/bin/bash
# 3-version benchmark on the public chr22 VCF (Test/1000genomes_chr22):
# PopLDdecay2 vs PopLDdecay-3.43 vs PopLDdecay-3.45. Compares wall/RSS per run
# and diff-checks .stat.gz content against the 3.45 reference (the version
# PopLDdecay2 was validated against).
# -T 16 on chr22 (single chromosome) exercises the fastpath position-shard path.
set -uo pipefail
cd "$(dirname "$0")/.."

BB=${BB:-Test/1000genomes_chr22/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz}
WORK=Test//ldbench3
mkdir -p "$WORK"
O43=/home/heweiming/01.Software/PopLDdecay-3.43/bin/PopLDdecay
O45=PopLDdecay-3.45/bin/PopLDdecay
NEW=bin/PopLDdecay2

wall_secs() { echo "$1" | awk -F: '{ if (NF==3) print $1*3600+$2*60+$3; else print $1*60+$2 }'; }

bench() {  # name bin extra...
  local name="$1" bin="$2"; shift 2
  /usr/bin/time -v "$bin" -InVCF "$BB" -OutStat "$WORK/$name.stat" "$@" \
    > "$WORK/$name.out" 2> "$WORK/$name.time"
  local rc=$?
  local w r
  w=$(wall_secs "$(grep 'Elapsed' "$WORK/$name.time" | awk '{print $8}')")
  r=$(grep 'Maximum resident' "$WORK/$name.time" | awk '{print $6}')
  echo "  $name: rc=$rc wall=${w}s rss=$((r/1024))MB"
}

echo "== version benchmark on $BB =="
echo "data: $(zcat "$BB" | grep -vc '^#') sites, $(zcat "$BB" | grep -m1 '^#CHROM' | cut -f10- | wc -w) samples"
date +%s > "$WORK/start.ts"

bench old.343  "$O43"
bench old.345  "$O45"
bench new.t1   "$NEW" -T 1
bench new.t8   "$NEW" -T 8
bench new.t16  "$NEW" -T 16

date +%s > "$WORK/end.ts"
echo "== result consistency (decompressed .stat.gz diff vs 3.45) =="
for name in old.343 new.t1 new.t8 new.t16; do
  d=$(diff <(gzip -cd "$WORK/old.345.stat.gz") <(gzip -cd "$WORK/$name.stat.gz") | wc -l)
  echo "  $name vs 3.45: diff=$d"
done
echo "== elapsed total: $(( $(cat "$WORK/end.ts") - $(cat "$WORK/start.ts") ))s =="
