#!/bin/bash
# Paper benchmark 2: chr22 1000 Genomes Phase 3 sub-population LD decay.
# Test/sample.group is the two-column `sample<TAB>superpop` file PopLDdecay2
# consumes natively (one run, five outputs). PopLDdecay-3.45 only takes one
# single-column list per run, so the five lists are split out and 3.45 is run
# once per population; the decompressed .stat.gz files must be byte-identical
# per population. Also measures wall/RSS of the single PopLDdecay2 run vs the
# five sequential 3.45 runs.
# NOTE: the chr22 VCF copied from the 1000 Genomes website ships with a STALE
# .tbi index (data rewritten after the index was built), which breaks the indexed
# fast path and makes PopLDdecay2 output differ from 3.45. Rebuild the index
# before running:
#   Test/1000genomes_chr22/*.vcf.gz.tbi   (delete the old .tbi)
#   /data/01.SoftWare/samtools-1.23/bin/tabix -f Test/1000genomes_chr22/*.vcf.gz
set -uo pipefail
cd "$(dirname "$0")/.."

VCF=${VCF:-Test/1000genomes_chr22/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz}
GRP=Test/sample.group
WORK=Test//paper_subpop
mkdir -p "$WORK"
OLD=PopLDdecay-3.45/bin/PopLDdecay
NEW=bin/PopLDdecay2
MAXD=${MAXD:-300}

wall_secs() { echo "$1" | awk -F: '{ if (NF==3) print $1*3600+$2*60+$3; else print $1*60+$2 }'; }

echo "== chr22 sub-population benchmark (MaxDist=${MAXD}kb) =="
echo "group: $GRP"

# split two-column group file into per-population single-column lists
for g in AFR AMR EAS EUR SAS; do
  awk -v g="$g" '$2==g {print $1}' "$GRP" > "$WORK/$g.list"
done
echo "sample counts: $(for g in AFR AMR EAS EUR SAS; do echo -n "$g=$(wc -l < "$WORK/$g.list") "; done)"

date +%s > "$WORK/start.ts"

# 3.45: one run per population (single-column -SubPop)
for g in AFR AMR EAS EUR SAS; do
  /usr/bin/time -v "$OLD" -InVCF "$VCF" -OutStat "$WORK/old.$g.stat" -MaxDist "$MAXD" -SubPop "$WORK/$g.list" \
    > "$WORK/old.$g.out" 2> "$WORK/old.$g.time"
  echo "  old.$g: rc=$? wall=$(wall_secs "$(grep 'Elapsed' "$WORK/old.$g.time" | awk '{print $8}')")s rss=$(($(grep 'Maximum resident' "$WORK/old.$g.time" | awk '{print $6}')/1024))MB"
done

# PopLDdecay2: single run, two-column file, 8 threads
/usr/bin/time -v "$NEW" -InVCF "$VCF" -OutStat "$WORK/new.stat" -MaxDist "$MAXD" -T 8 -SubPop "$GRP" \
  > "$WORK/new.out" 2> "$WORK/new.time"
echo "  new(t8): rc=$? wall=$(wall_secs "$(grep 'Elapsed' "$WORK/new.time" | awk '{print $8}')")s rss=$(($(grep 'Maximum resident' "$WORK/new.time" | awk '{print $6}')/1024))MB"

date +%s > "$WORK/end.ts"

echo "== per-population consistency (decompressed .stat.gz diff vs 3.45) =="
for g in AFR AMR EAS EUR SAS; do
  d=$(diff <(gzip -cd "$WORK/old.$g.stat.gz") <(gzip -cd "$WORK/new.$g.stat.gz") | wc -l)
  echo "  $g: diff=$d"
done
echo "== output files =="
ls -l "$WORK"/new.*.stat.gz "$WORK"/old.*.stat.gz | awk '{print "  "$NF" "$5}'
echo "== elapsed total: $(( $(cat "$WORK/end.ts") - $(cat "$WORK/start.ts") ))s =="
