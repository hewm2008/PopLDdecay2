#!/bin/bash
# Golden-diff: PopLDdecay2 output must be byte-identical (after gunzip) to the
# original PopLDdecay for the same input / parameters.
#
# Needs PopLDdecay 3.45 in the package root: PopLDdecay-3.45/bin/PopLDdecay
set -uo pipefail
cd "$(dirname "$0")/../.."

# Without a valid OLD every "== old" case compares two empty streams and passes,
# so refuse to start instead. Override with OLD=/path/to/PopLDdecay.
OLD=${OLD:-PopLDdecay-3.45/bin/PopLDdecay}
if [ ! -x "$OLD" ]; then
  echo "FATAL: reference PopLDdecay 3.45 not executable at '$OLD'"
  echo "       Set OLD=/path/to/PopLDdecay-3.45/bin/PopLDdecay and re-run."
  echo "       Refusing to continue: without it the byte-identity checks are vacuous."
  exit 2
fi
NEW=bin/PopLDdecay2
BGZIP=${BGZIP:-bgzip}
TABIX=${TABIX:-tabix}
WORK=${WORK:-$(mktemp -d)}
mkdir -p "$WORK"

fail=0

# cmp_gz LABEL OLD NEW -- decompressed content must be identical
cmp_gz() {
  # A missing file must not compare as two empty streams and pass: that is
  # the failure mode this script's header warns about.
  for f in "$2" "$3"; do
    if [ ! -f "$f" ]; then
      echo "FAIL  $1 (missing $f)"; fail=1; return
    fi
  done
  if diff <(gzip -cd "$2") <(gzip -cd "$3") > "$WORK/diff.last"; then
    echo "PASS  $1"
  else
    echo "FAIL  $1"; head -20 "$WORK/diff.last"; fail=1
  fi
}

# cmp_raw LABEL OLD NEW -- the files themselves must be identical
cmp_raw() {
  if diff "$2" "$3" > "$WORK/diff.last"; then
    echo "PASS  $1"
  else
    echo "FAIL  $1"; head -20 "$WORK/diff.last"; fail=1
  fi
}

# ---- inputs (plain / gzip / bgzip / phased) ----
python3 RunTest/06.Verify/gen_sim_vcf.py -o "$WORK/sim.vcf" --nsamples 60 --sites 2000 --chroms 4 >/dev/null
python3 RunTest/06.Verify/gen_sim_vcf.py -o "$WORK/sim.vcf.gz" --gzip --nsamples 60 --sites 2000 --chroms 4 >/dev/null
"$BGZIP" -f -c "$WORK/sim.vcf" > "$WORK/sim.bgzip.vcf.gz"
python3 RunTest/06.Verify/gen_sim_vcf.py -o "$WORK/sim_phased.vcf.gz" --gzip --phased --nsamples 60 --sites 2000 --chroms 4 >/dev/null

for name in sim.vcf sim.vcf.gz sim.bgzip.vcf.gz sim_phased.vcf.gz; do
  in="$WORK/$name"
  "$OLD" -InVCF "$in" -OutStat "$WORK/old.$name.stat" 2>/dev/null
  "$NEW" -InVCF "$in" -OutStat "$WORK/new.$name.stat" 2>/dev/null
  cmp_gz "$name" "$WORK/old.$name.stat.gz" "$WORK/new.$name.stat.gz"
done

# ---- parameter variation ----
"$OLD" -InVCF "$WORK/sim.vcf.gz" -OutStat "$WORK/old.param.stat" -MAF 0.1 -Het 0.5 -Miss 0.1 -MaxDist 100 2>/dev/null
"$NEW" -InVCF "$WORK/sim.vcf.gz" -OutStat "$WORK/new.param.stat" -MAF 0.1 -Het 0.5 -Miss 0.1 -MaxDist 100 2>/dev/null
cmp_gz "param(-MAF 0.1 -Het 0.5 -Miss 0.1 -MaxDist 100)" "$WORK/old.param.stat.gz" "$WORK/new.param.stat.gz"

# ---- default MaxDist (300kb) is covered by sim runs; also test -OutStat already ending in .stat.gz ----
"$OLD" -InVCF "$WORK/sim.vcf.gz" -OutStat "$WORK/old.prefix.stat.gz" 2>/dev/null
"$NEW" -InVCF "$WORK/sim.vcf.gz" -OutStat "$WORK/new.prefix.stat.gz" 2>/dev/null
cmp_gz "prefix(-OutStat *.stat.gz)" "$WORK/old.prefix.stat.gz" "$WORK/new.prefix.stat.gz"

# ---- multithreaded self-consistency: -T 8 must equal -T 1 (formatted) ----
"$NEW" -InVCF "$WORK/sim.vcf.gz" -OutStat "$WORK/t1.stat" -T 1 2>/dev/null
"$NEW" -InVCF "$WORK/sim.vcf.gz" -OutStat "$WORK/t8.stat" -T 8 2>/dev/null
cmp_gz "threads(-T 8 vs -T 1)" "$WORK/t1.stat.gz" "$WORK/t8.stat.gz"

# ---- indexed fast path: bgzip + tabix index, -T 4 must match old (fast path) ----
"$BGZIP" -f -c "$WORK/sim.vcf" > "$WORK/sim.idx.vcf.gz"
"$TABIX" -p vcf "$WORK/sim.idx.vcf.gz"
"$OLD" -InVCF "$WORK/sim.idx.vcf.gz" -OutStat "$WORK/old.simidx.stat" 2>/dev/null
"$NEW" -InVCF "$WORK/sim.idx.vcf.gz" -OutStat "$WORK/new.simidx.stat" -T 4 2>/dev/null
cmp_gz "fastpath indexed(bgzip+tabix, -T 4 vs old)" "$WORK/old.simidx.stat.gz" "$WORK/new.simidx.stat.gz"
# sharded fastpath: threads(8) > contigs(4) -> position shards; must still match old
"$OLD" -InVCF "$WORK/sim.idx.vcf.gz" -OutStat "$WORK/old.shard.stat" 2>/dev/null
"$NEW" -InVCF "$WORK/sim.idx.vcf.gz" -OutStat "$WORK/new.shard.stat" -T 8 2>/dev/null
cmp_gz "fastpath sharded(threads>contigs, -T 8 vs old)" "$WORK/old.shard.stat.gz" "$WORK/new.shard.stat.gz"
# phased indexed counterpart (for the fastpath .LD.gz checks below)
gzip -cd "$WORK/sim_phased.vcf.gz" | "$BGZIP" -f -c > "$WORK/sim_phased.idx.vcf.gz"
"$TABIX" -p vcf "$WORK/sim_phased.idx.vcf.gz"
# sharded fastpath .LD.gz byte-identical to the fallback (OutType 7, phased)
"$NEW" -InVCF "$WORK/sim_phased.idx.vcf.gz" -OutStat "$WORK/new.shardld.t1.stat" -OutType 7 -T 1 2>/dev/null
"$NEW" -InVCF "$WORK/sim_phased.idx.vcf.gz" -OutStat "$WORK/new.shardld.t8.stat" -OutType 7 -T 8 2>/dev/null
cmp_gz "fastpath sharded (OutType 7 .LD.gz T8 vs T1)" "$WORK/new.shardld.t1.LD.gz" "$WORK/new.shardld.t8.LD.gz"
# single-chromosome indexed fastpath (0.6.8): contigs>=1 now admitted; threads
# (8) > 1 chr -> position shards on the single chromosome; must match old
python3 RunTest/06.Verify/gen_sim_vcf.py -o "$WORK/sim1.vcf" --nsamples 60 --sites 1500 --chroms 1 >/dev/null
"$BGZIP" -f -c "$WORK/sim1.vcf" > "$WORK/sim1.idx.vcf.gz"
"$TABIX" -p vcf "$WORK/sim1.idx.vcf.gz"
"$OLD" -InVCF "$WORK/sim1.idx.vcf.gz" -OutStat "$WORK/old.sim1.stat" 2>/dev/null
"$NEW" -InVCF "$WORK/sim1.idx.vcf.gz" -OutStat "$WORK/new.sim1.stat" -T 8 2>/dev/null
cmp_gz "fastpath single-chr (bgzip+tabix, -T 8 vs old)" "$WORK/old.sim1.stat.gz" "$WORK/new.sim1.stat.gz"

# ---- OutType 0-8: .stat.gz (+ .LD.gz for 3/6/7/8) must match old, unphased & phased ----
for ot in 0 1 2 3 4 5 6 7 8; do
  for phase in unph phased; do
    if [ "$phase" = phased ]; then in="$WORK/sim_phased.vcf.gz"; else in="$WORK/sim.vcf.gz"; fi
    "$OLD" -InVCF "$in" -OutStat "$WORK/old.ot$ot.$phase.stat" -OutType "$ot" 2>/dev/null
    "$NEW" -InVCF "$in" -OutStat "$WORK/new.ot$ot.$phase.stat" -OutType "$ot" -T 1 2>/dev/null
    cmp_gz "OutType $ot $phase (.stat.gz)" "$WORK/old.ot$ot.$phase.stat.gz" "$WORK/new.ot$ot.$phase.stat.gz"
    case "$ot" in
      3|6|7|8)
        cmp_gz "OutType $ot $phase (.LD.gz)" "$WORK/old.ot$ot.$phase.LD.gz" "$WORK/new.ot$ot.$phase.LD.gz"
        # indexed fast path (-T 4) must produce the identical .LD.gz
        if [ "$phase" = phased ]; then idxin="$WORK/sim_phased.idx.vcf.gz"; else idxin="$WORK/sim.idx.vcf.gz"; fi
        "$NEW" -InVCF "$idxin" -OutStat "$WORK/new.ot$ot.fast.stat" -OutType "$ot" -T 4 2>/dev/null
        cmp_gz "OutType $ot $phase (fastpath .LD.gz)" "$WORK/new.ot$ot.$phase.LD.gz" "$WORK/new.ot$ot.fast.LD.gz"
        ;;
    esac
  done
done

# ---- -L region restrict (new feature, no old-binary target: self-consistency) ----
# (a) fallback vs fastpath on the indexed input must be formatted-identical
"$NEW" -InVCF "$WORK/sim.idx.vcf.gz" -OutStat "$WORK/l1.t1.stat" -L chr1 -T 1 2>/dev/null
"$NEW" -InVCF "$WORK/sim.idx.vcf.gz" -OutStat "$WORK/l1.t4.stat" -L chr1 -T 4 2>/dev/null
cmp_gz "-L chr1 (fallback T1 vs fastpath T4)" "$WORK/l1.t1.stat.gz" "$WORK/l1.t4.stat.gz"
# (b) -L chr1 on the full VCF == running on a chr1-only sub-VCF (byte-identical)
awk -F'\t' '$1=="#CHROM" || $1=="chr1"' "$WORK/sim.vcf" > "$WORK/sim.chr1.vcf"
"$NEW" -InVCF "$WORK/sim.vcf" -OutStat "$WORK/eq.full.stat" -L chr1 -T 1 2>/dev/null
"$NEW" -InVCF "$WORK/sim.chr1.vcf" -OutStat "$WORK/eq.sub.stat" -T 1 2>/dev/null
cmp_gz "-L chr1 == chr1-only sub-VCF" "$WORK/eq.full.stat.gz" "$WORK/eq.sub.stat.gz"
# (c) BED file (0-based half-open) fallback vs fastpath
printf 'chr1\t999\t1000000\nchr2\t1\t5000000\n' > "$WORK/regions.bed"
"$NEW" -InVCF "$WORK/sim.idx.vcf.gz" -OutStat "$WORK/bed.t1.stat" -L "$WORK/regions.bed" -T 1 2>/dev/null
"$NEW" -InVCF "$WORK/sim.idx.vcf.gz" -OutStat "$WORK/bed.t4.stat" -L "$WORK/regions.bed" -T 4 2>/dev/null
cmp_gz "-L BED (fallback T1 vs fastpath T4)" "$WORK/bed.t1.stat.gz" "$WORK/bed.t4.stat.gz"
# (d) empty region (chr not present) -> header-only output, both paths equal
"$NEW" -InVCF "$WORK/sim.idx.vcf.gz" -OutStat "$WORK/empty.t1.stat" -L chrY -T 1 2>/dev/null
"$NEW" -InVCF "$WORK/sim.idx.vcf.gz" -OutStat "$WORK/empty.t4.stat" -L chrY -T 4 2>/dev/null
cmp_gz "-L empty region (header only)" "$WORK/empty.t1.stat.gz" "$WORK/empty.t4.stat.gz"

# ---- B5: indexed VCF with data on a contig that has NO ##contig line ----
# The fastpath task list comes from the tabix index (tbx_seqnames); a
# header-only list would silently drop chr2's pairs (the original single-pass
# reader sees every contig). All three runs must agree.
grep -v '^##contig=<ID=chr2,' "$WORK/sim.vcf" > "$WORK/undec.vcf"
"$BGZIP" -f -c "$WORK/undec.vcf" > "$WORK/undec.vcf.gz"
"$TABIX" -p vcf "$WORK/undec.vcf.gz"
"$NEW" -InVCF "$WORK/undec.vcf.gz" -OutStat "$WORK/undec.t1.stat" -T 1 2>/dev/null
"$NEW" -InVCF "$WORK/undec.vcf.gz" -OutStat "$WORK/undec.t4.stat" -T 4 2>/dev/null
"$OLD" -InVCF "$WORK/undec.vcf.gz" -OutStat "$WORK/undec.old.stat" 2>/dev/null
if diff <(gzip -cd "$WORK/undec.t1.stat.gz") <(gzip -cd "$WORK/undec.t4.stat.gz") > "$WORK/diff.undec.t" &&
   diff <(gzip -cd "$WORK/undec.t4.stat.gz") <(gzip -cd "$WORK/undec.old.stat.gz") > "$WORK/diff.undec.old"; then
  echo "PASS  undeclared ##contig (T1 vs T4 vs old)"
else
  echo "FAIL  undeclared ##contig"; head -20 "$WORK/diff.undec.t"; head -20 "$WORK/diff.undec.old"; fail=1
fi

# ---- B6: duplicate -L region must not double-count on the fast path ----
# The fallback folds identical regions into one label-keyed buffer; without
# dedup the fastpath would run one task per duplicate (pairs doubled).
"$NEW" -InVCF "$WORK/sim.idx.vcf.gz" -OutStat "$WORK/dupreg.t1.stat" -L chr1,chr1 -T 1 2>/dev/null
"$NEW" -InVCF "$WORK/sim.idx.vcf.gz" -OutStat "$WORK/dupreg.t4.stat" -L chr1,chr1 -T 4 2>/dev/null
cmp_gz "duplicate -L chr1,chr1 (T1 vs T4)" "$WORK/dupreg.t1.stat.gz" "$WORK/dupreg.t4.stat.gz"

# ---- B8: ':'-containing contig name -> fastpath declines (tabix query unsafe),
#      single-pass fallback handles it; T1 vs T4 vs old must agree ----
sed 's/\bchr2\b/chr2:alt/g' "$WORK/sim.vcf" > "$WORK/colon.vcf"
"$BGZIP" -f -c "$WORK/colon.vcf" > "$WORK/colon.vcf.gz"
"$TABIX" -p vcf "$WORK/colon.vcf.gz"
"$NEW" -InVCF "$WORK/colon.vcf.gz" -OutStat "$WORK/colon.t1.stat" -T 1 2>/dev/null
"$NEW" -InVCF "$WORK/colon.vcf.gz" -OutStat "$WORK/colon.t4.stat" -T 4 2>/dev/null
"$OLD" -InVCF "$WORK/colon.vcf.gz" -OutStat "$WORK/colon.old.stat" 2>/dev/null
if diff <(gzip -cd "$WORK/colon.t1.stat.gz") <(gzip -cd "$WORK/colon.t4.stat.gz") > "$WORK/diff.colon.t" &&
   diff <(gzip -cd "$WORK/colon.t4.stat.gz") <(gzip -cd "$WORK/colon.old.stat.gz") > "$WORK/diff.colon.old"; then
  echo "PASS  ':'-contig fastpath decline (T1 vs T4 vs old)"
else
  echo "FAIL  ':'-contig"; head -20 "$WORK/diff.colon.t"; head -20 "$WORK/diff.colon.old"; fail=1
fi

# ---- B7: single-char GTs in the MIDDLE sample columns (last column kept 3-char,
#      so the original is fully deterministic: split3 raw reads go into the next
#      field for short tokens) -> byte-identical to the original ----
python3 - "$WORK/sim.vcf" "$WORK/haploid_mid.vcf" <<'PYEOF'
import sys
inp, out = sys.argv[1], sys.argv[2]
with open(inp) as f, open(out, 'w') as g:
    for line in f:
        if line.startswith('#'):
            g.write(line); continue
        p = line.rstrip('\n').split('\t')
        for i in range(9, len(p) - 1):
            if p[i] in ('0/0', '1/1', '0/1', '1/0'):
                p[i] = p[i][0]
        g.write('\t'.join(p) + '\n')
PYEOF
"$OLD" -InVCF "$WORK/haploid_mid.vcf" -OutStat "$WORK/haploid_mid.old.stat" 2>/dev/null
"$NEW" -InVCF "$WORK/haploid_mid.vcf" -OutStat "$WORK/haploid_mid.new.stat" -T 1 2>/dev/null
cmp_gz "single-char GT (mid columns, vs old)" "$WORK/haploid_mid.old.stat.gz" "$WORK/haploid_mid.new.stat.gz"

# ---- -SubPop subset (whitespace-separated list, 30 of 60 samples) ----
{ for i in $(seq 1 2 59); do echo -n "S$i "; done; echo; echo "# comment"; echo "S99"; } > "$WORK/pop1.list"
for phase in unphased phased; do
  if [ "$phase" = phased ]; then in="$WORK/sim_phased.vcf.gz"; else in="$WORK/sim.vcf.gz"; fi
  "$OLD" -InVCF "$in" -OutStat "$WORK/old.sub.$phase.stat" -SubPop "$WORK/pop1.list" 2>/dev/null
  "$NEW" -InVCF "$in" -OutStat "$WORK/new.sub.$phase.stat" -SubPop "$WORK/pop1.list" -T 1 2>/dev/null
  cmp_gz "-SubPop $phase (30/60 vs old)" "$WORK/old.sub.$phase.stat.gz" "$WORK/new.sub.$phase.stat.gz"
done
# SubPop on the indexed input: fastpath (T4) must equal fallback (T1)
"$NEW" -InVCF "$WORK/sim.idx.vcf.gz" -OutStat "$WORK/sub.t1.stat" -SubPop "$WORK/pop1.list" -T 1 2>/dev/null
"$NEW" -InVCF "$WORK/sim.idx.vcf.gz" -OutStat "$WORK/sub.t4.stat" -SubPop "$WORK/pop1.list" -T 4 2>/dev/null
cmp_gz "-SubPop indexed (fallback T1 vs fastpath T4)" "$WORK/sub.t1.stat.gz" "$WORK/sub.t4.stat.gz"
# ---- two-column -SubPop (sample<TAB>groupid): each group must equal a
#      single-column run of that subset (self-consistency; the original tool has
#      no two-column mode) ----
awk '{ for (i=1;i<=NF;i++) print $i "\t" (i%2==1?"popA":"popB") }' "$WORK/pop1.list" > "$WORK/twocol.list"
awk '$2=="popA"{print $1}' "$WORK/twocol.list" > "$WORK/popA.list"
awk '$2=="popB"{print $1}' "$WORK/twocol.list" > "$WORK/popB.list"
"$NEW" -InVCF "$WORK/sim.vcf.gz" -OutStat "$WORK/tc.stat" -SubPop "$WORK/twocol.list" -T 1 2>/dev/null
"$NEW" -InVCF "$WORK/sim.vcf.gz" -OutStat "$WORK/tcA.stat" -SubPop "$WORK/popA.list" -T 1 2>/dev/null
"$NEW" -InVCF "$WORK/sim.vcf.gz" -OutStat "$WORK/tcB.stat" -SubPop "$WORK/popB.list" -T 1 2>/dev/null
ok=1
for f in "$WORK/tc.popA.stat.gz" "$WORK/tc.popB.stat.gz"; do
  [ -f "$f" ] || { echo "FAIL  two-column -SubPop (missing $f)"; ok=0; fail=1; }
done
# A single single-column file is legacy mode -> {prefix}.stat.gz (no label).
if [ "$ok" = 1 ]; then
  if diff <(gzip -cd "$WORK/tc.popA.stat.gz") <(gzip -cd "$WORK/tcA.stat.gz") > "$WORK/diff.tcA" \
     && diff <(gzip -cd "$WORK/tc.popB.stat.gz") <(gzip -cd "$WORK/tcB.stat.gz") > "$WORK/diff.tcB"; then
    echo "PASS  two-column -SubPop (per-group == single-column run)"
  else
    echo "FAIL  two-column -SubPop"; head -10 "$WORK/diff.tcA" "$WORK/diff.tcB"; fail=1
  fi
fi
# repeated -SubPop (one single-column file per group) equals the two-column run
"$NEW" -InVCF "$WORK/sim.vcf.gz" -OutStat "$WORK/rep.stat" -SubPop "$WORK/popA.list" -SubPop "$WORK/popB.list" -T 1 2>/dev/null
if diff <(gzip -cd "$WORK/rep.popA.stat.gz") <(gzip -cd "$WORK/tc.popA.stat.gz") > "$WORK/diff.repA" \
   && diff <(gzip -cd "$WORK/rep.popB.stat.gz") <(gzip -cd "$WORK/tc.popB.stat.gz") > "$WORK/diff.repB"; then
  echo "PASS  repeated -SubPop (per-group == single-column run)"
else
  echo "FAIL  repeated -SubPop"; head -10 "$WORK/diff.repA" "$WORK/diff.repB"; fail=1
fi

# ---- -OutFilterSNP (plain-gzip filter file, decompressed == old, raw bytes too) ----
"$OLD" -InVCF "$WORK/sim.vcf.gz" -OutStat "$WORK/old.flt.stat" -OutFilterSNP 2>/dev/null
"$NEW" -InVCF "$WORK/sim.vcf.gz" -OutStat "$WORK/new.flt.stat" -OutFilterSNP 2>/dev/null
ok=1
for f in old.flt new.flt; do
  [ -f "$WORK/$f.vcf.filter.gz" ] || { echo "FAIL  -OutFilterSNP ($f missing)"; ok=0; fail=1; }
done
if [ "$ok" = 1 ]; then
  if diff <(gzip -cd "$WORK/old.flt.vcf.filter.gz") <(gzip -cd "$WORK/new.flt.vcf.filter.gz") > "$WORK/diff.flt" \
     && cmp -s "$WORK/old.flt.vcf.filter.gz" "$WORK/new.flt.vcf.filter.gz"; then
    echo "PASS  -OutFilterSNP (filter.gz == old, raw bytes too)"
  else
    echo "FAIL  -OutFilterSNP"; head -20 "$WORK/diff.flt"; fail=1
  fi
fi
# SubPop + OutFilterSNP writes NO filter file (as the original)
"$NEW" -InVCF "$WORK/sim.vcf.gz" -OutStat "$WORK/subflt.stat" -SubPop "$WORK/pop1.list" -OutFilterSNP 2>/dev/null
if [ ! -f "$WORK/subflt.vcf.filter.gz" ]; then
  echo "PASS  -SubPop+-OutFilterSNP (no filter file, matches original)"
else
  echo "FAIL  -SubPop+-OutFilterSNP (unexpected filter file)"; fail=1
fi

# ---- real 1000G chr22 head accuracy (optional: skipped if data absent) ----
# Uses the first 2000 sites of the public 1000 Genomes chr22 VCF (2,504 samples,
# phased) - fast enough for the per-build gate while exercising the big-N path.
# 1000 Genomes Phase 3 chr22. Point CHR22_SRC at a local copy to enable the
# real-data cases below (they are skipped when the file is absent); the shipped
# .tbi may be stale after a re-download, so the derived subsets are re-indexed.
CHR22_SRC=${CHR22_SRC:-RunTest/00.GetData/chr22/chr22.vcf.gz}
if [ -f "$CHR22_SRC" ]; then
  if [ ! -f "$WORK/chr22head.vcf.gz" ]; then
    gzip -cd "$CHR22_SRC" | awk 'BEGIN{k=0} /^#/{print; next} {k++; if(k<=2000) print; else exit}' > "$WORK/chr22head.vcf"
    "$BGZIP" -f -c "$WORK/chr22head.vcf" > "$WORK/chr22head.vcf.gz" && rm -f "$WORK/chr22head.vcf"
    "$TABIX" -p vcf "$WORK/chr22head.vcf.gz"
  fi
  for T in 1 4 8; do
    "$OLD" -InVCF "$WORK/chr22head.vcf.gz" -OutStat "$WORK/old.c22.t$T.stat" 2>/dev/null
    "$NEW" -InVCF "$WORK/chr22head.vcf.gz" -OutStat "$WORK/new.c22.t$T.stat" -T "$T" 2>/dev/null
    cmp_gz "real chr22 head -T $T" "$WORK/old.c22.t$T.stat.gz" "$WORK/new.c22.t$T.stat.gz"
  done
  # SubPop: first 500 of the 2504 samples (matches old; fastpath == fallback)
  gunzip -c "$CHR22_SRC" | grep -m1 "^#CHROM" | cut -f10- | tr '\t' '\n' \
    | head -500 | tr '\n' ' ' > "$WORK/c22.g1.list"
  "$OLD" -InVCF "$WORK/chr22head.vcf.gz" -OutStat "$WORK/old.c22.sub.stat" -SubPop "$WORK/c22.g1.list" 2>/dev/null
  "$NEW" -InVCF "$WORK/chr22head.vcf.gz" -OutStat "$WORK/new.c22.sub.stat" -SubPop "$WORK/c22.g1.list" -T 1 2>/dev/null
  cmp_gz "real chr22 head -SubPop (500 samples vs old)" "$WORK/old.c22.sub.stat.gz" "$WORK/new.c22.sub.stat.gz"
  "$NEW" -InVCF "$WORK/chr22head.vcf.gz" -OutStat "$WORK/new.c22.subt.stat" -SubPop "$WORK/c22.g1.list" -T 8 2>/dev/null
  cmp_gz "real chr22 head -SubPop indexed (fallback T1 vs fastpath T8)" "$WORK/new.c22.sub.stat.gz" "$WORK/new.c22.subt.stat.gz"
else
  echo "SKIP  real chr22 head ($CHR22_SRC not present)"
fi

# ---- -EHH (M5): raw .ehh.gz bytes must equal the original ogzstream output ----
ehh_ok() {  # name in core [extra...]
  local name="$1" in="$2" core="$3"; shift 3
  "$OLD" -InVCF "$in" -OutStat "$WORK/old.$name.stat" -EHH "$core" "$@" >/dev/null 2>&1
  "$NEW" -InVCF "$in" -OutStat "$WORK/new.$name.stat" -EHH "$core" "$@" >/dev/null 2>&1
  if cmp -s "$WORK/old.$name.ehh.gz" "$WORK/new.$name.ehh.gz" \
     && diff <(gzip -cd "$WORK/old.$name.ehh.gz") <(gzip -cd "$WORK/new.$name.ehh.gz") >/dev/null; then
    echo "PASS  -EHH $name (raw .ehh.gz bytes == old)"
  else
    echo "FAIL  -EHH $name (raw bytes differ)"; fail=1
  fi
  rm -f "$WORK/old.$name.stat.gz" "$WORK/new.$name.stat.gz" \
       "$WORK/old.$name.ehh.gz" "$WORK/new.$name.ehh.gz" \
       "$WORK/old.$name.ehh.pdf" "$WORK/new.$name.ehh.pdf" \
       "$WORK/old.$name.ehh.png" "$WORK/new.$name.ehh.png"
}
# Error paths leave an empty 20-byte .ehh.gz (missing chr / <8 no-missing SNPs).
ehh_ok ehh.nochr "$WORK/sim.vcf.gz" chrXX:700000
ehh_ok ehh.few   "$WORK/sim.vcf.gz" chr1:700000
if [ -f "$CHR22_SRC" ]; then
  [ -f "$WORK/c22.g1.list" ] || { gunzip -c "$CHR22_SRC" | grep -m1 "^#CHROM" \
    | cut -f10- | tr '\t' '\n' | head -500 | tr '\n' ' ' > "$WORK/c22.g1.list"; }
  ehh_ok ehh.c22     "$WORK/chr22head.vcf.gz" 22:$(gzip -cd "$WORK/chr22head.vcf.gz" | grep -v '^#' | awk 'NR==1000{print $2}')
  ehh_ok ehh.c22.ot8 "$WORK/chr22head.vcf.gz" 22:$(gzip -cd "$WORK/chr22head.vcf.gz" | grep -v '^#' | awk 'NR==1000{print $2}') -OutType 8
  ehh_ok ehh.c22.sub "$WORK/chr22head.vcf.gz" 22:$(gzip -cd "$WORK/chr22head.vcf.gz" | grep -v '^#' | awk 'NR==1000{print $2}') -SubPop "$WORK/c22.g1.list"
else
  echo "SKIP  -EHH real chr22 head (Test data not present)"
fi

# ---- contig with a very long name -----------------------------------------
# The .LD.gz and .ehh.gz rows start with the chromosome name. The original
# writes them through a C++ stream, which has no length limit; this tool used
# to format them into a fixed stack buffer and hand the length snprintf()
# returned (the length the row WOULD have had) to gzwrite(), so the tail of
# every row was read from beyond the buffer. With a long name the two must
# agree byte for byte.
LONG_NAME="chr$(printf 'L%.0s' {1..197})"
python3 RunTest/06.Verify/gen_sim_vcf.py -o "$WORK/longname.vcf" \
        --chrom-name "$LONG_NAME" --nsamples 60 --sites 400 --chroms 1 >/dev/null
for ot in 3 8; do
  "$OLD" -InVCF "$WORK/longname.vcf" -OutStat "$WORK/old.ln$ot" -OutType "$ot" -MaxDist 100 2>/dev/null
  "$NEW" -InVCF "$WORK/longname.vcf" -OutStat "$WORK/new.ln$ot" -OutType "$ot" -MaxDist 100 2>/dev/null
  cmp_gz "long contig name, OutType $ot (.LD.gz == old)" \
         "$WORK/old.ln$ot.LD.gz" "$WORK/new.ln$ot.LD.gz"
done
# the same input on the indexed fast path
"$BGZIP" -f -c "$WORK/longname.vcf" > "$WORK/longname.idx.vcf.gz"
"$TABIX" -p vcf "$WORK/longname.idx.vcf.gz"
"$OLD" -InVCF "$WORK/longname.idx.vcf.gz" -OutStat "$WORK/old.lnfp" -OutType 8 -MaxDist 100 2>/dev/null
"$NEW" -InVCF "$WORK/longname.idx.vcf.gz" -OutStat "$WORK/new.lnfp" -OutType 8 -MaxDist 100 -T 4 2>/dev/null
cmp_gz "long contig name, indexed fast path (.LD.gz == old)" \
       "$WORK/old.lnfp.LD.gz" "$WORK/new.lnfp.LD.gz"
# EHH: the rows of the .ehh.gz table start with the chromosome name too
python3 RunTest/06.Verify/gen_sim_vcf.py -o "$WORK/longname_ehh.vcf" \
        --chrom-name "$LONG_NAME" --nsamples 60 --sites 20000 --chroms 1 >/dev/null
# The generator emits missing calls; -EHH drops every site with a missing
# sample and needs more than 8 fully called sites in the window.
awk 'BEGIN{FS=OFS="\t"} /^#/{print; next}
     { for (i=10; i<=NF; i++) gsub(/\./, "0", $i); print }' "$WORK/longname_ehh.vcf" \
  | "$BGZIP" -c > "$WORK/longname_ehh.vcf.gz"
"$TABIX" -p vcf "$WORK/longname_ehh.vcf.gz"
core=$(zcat "$WORK/longname_ehh.vcf.gz" | grep -v '^#' | awk 'NR==200{print $2}')
if [ -n "$core" ]; then
  "$OLD" -InVCF "$WORK/longname_ehh.vcf.gz" -OutStat "$WORK/old.lnehh" -EHH "$LONG_NAME:$core" >/dev/null 2>&1
  "$NEW" -InVCF "$WORK/longname_ehh.vcf.gz" -OutStat "$WORK/new.lnehh" -EHH "$LONG_NAME:$core" >/dev/null 2>&1
  if cmp -s "$WORK/old.lnehh.ehh.gz" "$WORK/new.lnehh.ehh.gz" \
     && diff <(gzip -cd "$WORK/old.lnehh.ehh.gz") \
             <(gzip -cd "$WORK/new.lnehh.ehh.gz") >/dev/null; then
    echo "PASS  long contig name, -EHH (raw .ehh.gz bytes == old)"
  else
    echo "FAIL  long contig name, -EHH (raw bytes differ)"; fail=1
  fi
  rm -f "$WORK"/old.lnehh.* "$WORK"/new.lnehh.*
fi

# ---- -InGenotype (M6): stat content must match old; filter.gz raw bytes ----
geno_ok() {  # name infile extra...
  local name="$1" in="$2"; shift 2
  "$OLD" -InGenotype "$in" -OutStat "$WORK/old.$name.stat" "$@" >/dev/null 2>&1
  "$NEW" -InGenotype "$in" -OutStat "$WORK/new.$name.stat" "$@" >/dev/null 2>&1
  cmp_gz "-InGenotype $name (.stat.gz == old)" "$WORK/old.$name.stat.gz" "$WORK/new.$name.stat.gz"
}
python3 RunTest/06.Verify/gen_sim_genotype.py -i "$WORK/sim.vcf.gz" -o "$WORK/sim.geno" --with-header >/dev/null
geno_ok geno "$WORK/sim.geno"
# gz input (plain-gzip genotype input)
gzip -c "$WORK/sim.geno" > "$WORK/sim.geno.gz"
geno_ok geno.gz "$WORK/sim.geno.gz"
# ---- bin/mis/plink2genotype.pl --------------------------------------------
# The PED/MAP converter had no test at all: everything above feeds -InGenotype
# from gen_sim_genotype.py, which bypasses the converter, so it could have
# broken without a single case noticing. gen_sim_ped.py writes the PED/MAP and
# the .geno they must convert to, giving an oracle instead of a comparison of
# the converter against itself.
#
# The converter is also run under a memory ceiling. It used to hold the whole
# PED as one Perl scalar per field, i.e. nsamples x (6 + 2*nvariants) scalars,
# which reached 15 GiB on a 1000 x 100000 fixture; the rewrite streams it.
p2g_one() {   # name prefix [converter args...]
  local name=$1 pre=$2; shift 2
  perl bin/mis/plink2genotype.pl -inPED "$pre.ped" -inMAP "$pre.map" \
       -outGenotype "$WORK/p2g.$name" "$@" >/dev/null 2>&1
  if [ ! -s "$WORK/p2g.$name" ]; then
    echo "FAIL  plink2genotype $name (no output)"; fail=1; return
  fi
  if cmp -s "$WORK/p2g.$name" "$pre.geno"; then
    echo "PASS  plink2genotype $name (== oracle .geno)"
  else
    echo "FAIL  plink2genotype $name"; diff "$pre.geno" "$WORK/p2g.$name" | head -5; fail=1
  fi
}
for p2genc in numeric letter; do
  python3 RunTest/06.Verify/gen_sim_ped.py -i "$WORK/sim.vcf.gz" -o "$WORK/p2g_$p2genc" \
          --encoding "$p2genc" >/dev/null 2>&1
  p2g_one "$p2genc" "$WORK/p2g_$p2genc"
  # The converter already reproduces the oracle byte for byte, so feeding both
  # to -InGenotype would compare a file with itself. What is worth checking is
  # that the converted file is a valid input that PopLDdecay 3.45 and
  # PopLDdecay2 read identically.
  "$OLD" -InGenotype "$WORK/p2g.$p2genc" -OutStat "$WORK/p2g.$p2genc.old" >/dev/null 2>&1
  "$NEW" -InGenotype "$WORK/p2g.$p2genc" -OutStat "$WORK/p2g.$p2genc.new" >/dev/null 2>&1
  if [ -s "$WORK/p2g.$p2genc.old.stat.gz" ] && \
     diff <(gzip -cd "$WORK/p2g.$p2genc.old.stat.gz") \
          <(gzip -cd "$WORK/p2g.$p2genc.new.stat.gz") > "$WORK/diff.p2g.$p2genc"; then
    echo "PASS  -InGenotype on converted PED ($p2genc, == 3.45)"
  else
    echo "FAIL  -InGenotype on converted PED ($p2genc)"; head -5 "$WORK/diff.p2g.$p2genc"; fail=1
  fi
  # Blocking must not change a single byte. The block is a sensible fraction of
  # the variant count: because the PED is re-read once per block, the cost grows
  # with variants/block, and a deliberately tiny block (7) took 390 s here and
  # turned the whole suite from minutes into a quarter of an hour.
  perl bin/mis/plink2genotype.pl -inPED "$WORK/p2g_$p2genc.ped" -inMAP "$WORK/p2g_$p2genc.map" \
       -outGenotype "$WORK/p2g.$p2genc.blk" -block 1000 >/dev/null 2>&1
  if [ -s "$WORK/p2g.$p2genc.blk" ] && cmp -s "$WORK/p2g.$p2genc.blk" "$WORK/p2g_$p2genc.geno"; then
    echo "PASS  plink2genotype -block 1000 ($p2genc, == single pass)"
  else
    echo "FAIL  plink2genotype -block 1000 ($p2genc)"; fail=1
  fi
done
# hand-built edge fixture: every merge_alleles branch, including allele pairs
# that have no IUPAC code (the converter concatenates those)
python3 RunTest/06.Verify/gen_sim_ped.py --edge -o "$WORK/p2g_edge" --encoding letter >/dev/null 2>&1
p2g_one edge "$WORK/p2g_edge"
# -maxMem is a ceiling, not a target: the requested ceiling must be respected
if [ -s "$WORK/p2g_numeric.ped" ]; then
  /usr/bin/time -f "%M" -o "$WORK/p2g.rss" \
    perl bin/mis/plink2genotype.pl -inPED "$WORK/p2g_numeric.ped" \
         -inMAP "$WORK/p2g_numeric.map" -outGenotype "$WORK/p2g.cap" \
         -maxMem 64 >/dev/null 2>&1
  rss=$(cat "$WORK/p2g.rss" 2>/dev/null)
  if [ -n "$rss" ] && [ "$rss" -le $((64 * 1024 * 110 / 100)) ]; then
    echo "PASS  plink2genotype -maxMem 64 (peak RSS $((rss / 1024)) MiB, within the ceiling)"
  else
    echo "FAIL  plink2genotype -maxMem 64 (peak RSS ${rss:-?} KB)"; fail=1
  fi
fi

# -OutFilterSNP: .genotype.filter.gz raw bytes == old (plain gzip) + stat
"$OLD" -InGenotype "$WORK/sim.geno" -OutStat "$WORK/old.genoflt.stat" -OutFilterSNP >/dev/null 2>&1
"$NEW" -InGenotype "$WORK/sim.geno" -OutStat "$WORK/new.genoflt.stat" -OutFilterSNP >/dev/null 2>&1
if [ -f "$WORK/old.genoflt.genotype.filter.gz" ] && [ -f "$WORK/new.genoflt.genotype.filter.gz" ] \
   && diff <(gzip -cd "$WORK/old.genoflt.genotype.filter.gz") <(gzip -cd "$WORK/new.genoflt.genotype.filter.gz") > "$WORK/diff.genoflt" \
   && cmp -s "$WORK/old.genoflt.genotype.filter.gz" "$WORK/new.genoflt.genotype.filter.gz" \
   && diff <(gzip -cd "$WORK/old.genoflt.stat.gz") <(gzip -cd "$WORK/new.genoflt.stat.gz") >/dev/null; then
  echo "PASS  -InGenotype -OutFilterSNP (filter.gz raw bytes == old)"
else
  echo "FAIL  -InGenotype -OutFilterSNP"; fail=1
fi
# SubPop (odd samples; header sample order == VCF order)
head -1 "$WORK/sim.geno" | cut -f3- | tr ' ' '\n' | awk 'NR%2==1' | tr '\n' ' ' > "$WORK/geno.g1.list"
geno_ok geno.sub "$WORK/sim.geno" -SubPop "$WORK/geno.g1.list"
# SubPop + singleton / tri-allelic sites -> BadSite warning, stat == old
printf '#CHROM\tPOS\tS1 S2 S3 S4\nchr1\t100\tA A A A\nchr1\t200\tA C G A\nchr1\t300\tA C A C\n' > "$WORK/tri.geno"
printf 'S1 S2 S3 S4\n' > "$WORK/tri.list"
geno_ok geno.tri "$WORK/tri.geno" -SubPop "$WORK/tri.list"
# VCF + genotype combined (merged SNPList) vs old
"$OLD" -InVCF "$WORK/sim.vcf.gz" -InGenotype "$WORK/sim.geno" -OutStat "$WORK/old.genomix.stat" >/dev/null 2>&1
"$NEW" -InVCF "$WORK/sim.vcf.gz" -InGenotype "$WORK/sim.geno" -OutStat "$WORK/new.genomix.stat" >/dev/null 2>&1
cmp_gz "-InGenotype mixed with -InVCF (.stat.gz == old)" "$WORK/old.genomix.stat.gz" "$WORK/new.genomix.stat.gz"
# EHH + genotype: raw .ehh.gz bytes == old
"$OLD" -InGenotype "$WORK/sim.geno" -OutStat "$WORK/old.genoehh.stat" -EHH chr1:2310 >/dev/null 2>&1
"$NEW" -InGenotype "$WORK/sim.geno" -OutStat "$WORK/new.genoehh.stat" -EHH chr1:2310 >/dev/null 2>&1
if cmp -s "$WORK/old.genoehh.ehh.gz" "$WORK/new.genoehh.ehh.gz" \
   && diff <(gzip -cd "$WORK/old.genoehh.ehh.gz") <(gzip -cd "$WORK/new.genoehh.ehh.gz") >/dev/null; then
  echo "PASS  -InGenotype -EHH (raw .ehh.gz bytes == old)"
else
  echo "FAIL  -InGenotype -EHH"; fail=1
fi
rm -f "$WORK/old.genoehh.stat.gz" "$WORK/new.genoehh.stat.gz" \
     "$WORK/old.genoehh.ehh.gz" "$WORK/new.genoehh.ehh.gz"

# ---- -InPlink (M7): native .bed/.bim/.fam must equal the same data as -InVCF ----
# PopLDdecay 3.45 has no PLINK input, so there is no reference binary here: the
# oracle is the identical dataset in VCF form, which IS covered against 3.45
# above. gen_sim_plink.py writes both representations from one allele matrix.
python3 RunTest/06.Verify/gen_sim_plink.py -o "$WORK/plk" --nsamples 61 --sites 4000 \
        --chroms 4 --seed 11 --with-header >/dev/null
plink_ok() {  # name infile extra...
  local name="$1"; shift
  "$NEW" -InPlink "$WORK/plk.bed" -OutStat "$WORK/plk.$name" "$@" >/dev/null 2>&1
  "$NEW" -InVCF "$WORK/plk.vcf" -OutStat "$WORK/vcf.$name" "$@" >/dev/null 2>&1
  cmp_gz "-InPlink $name (== -InVCF)" "$WORK/plk.$name.stat.gz" "$WORK/vcf.$name.stat.gz"
}
plink_ok default
plink_ok params -MAF 0.1 -Het 0.5 -Miss 0.1 -MaxDist 100
plink_ok outtype2 -OutType 2
plink_ok outtype8 -OutType 8
plink_ok threads -T 8
# OutType 8 also writes the pairwise .LD.gz
cmp_gz "-InPlink outtype8 .LD.gz (== -InVCF)" "$WORK/plk.outtype8.LD.gz" "$WORK/vcf.outtype8.LD.gz"
# -T 1 vs -T 8 on the PLINK input itself (the .bed is decoded in parallel chunks)
"$NEW" -InPlink "$WORK/plk.bed" -OutStat "$WORK/plk.t1" -T 1 >/dev/null 2>&1
cmp_gz "-InPlink -T 8 (parallel .bed decode == -T 1)" "$WORK/plk.t1.stat.gz" "$WORK/plk.threads.stat.gz"
# -SubPop: single-column (group per file) and two-column (group per distinct id)
awk 'NR>2 {print $2}' "$WORK/plk.fam" | awk 'NR%2==1' | tr '\n' ' ' > "$WORK/plk.g1.list"
plink_ok subpop -SubPop "$WORK/plk.g1.list"
awk 'NR>2 {printf "%s\t%s\n", $2, ((NR%3)==1 ? "A" : ((NR%3)==2 ? "B" : "C"))}' \
    "$WORK/plk.fam" | awk '!($2=="C")' > "$WORK/plk.grp.list"
for g in A B; do
  "$NEW" -InPlink "$WORK/plk.bed" -OutStat "$WORK/plk.g$g" -SubPop "$WORK/plk.grp.list" >/dev/null 2>&1
  "$NEW" -InVCF "$WORK/plk.vcf" -OutStat "$WORK/vcf.g$g" -SubPop "$WORK/plk.grp.list" >/dev/null 2>&1
  cmp_gz "-InPlink -SubPop group $g (== -InVCF)" "$WORK/plk.g$g.$g.stat.gz" "$WORK/vcf.g$g.$g.stat.gz"
done
# -L regions (per-region buffers must exclude cross-region pairs on both paths)
plink_ok region -L 1:2000-200000,2
# -EHH: raw .ehh.gz bytes
"$NEW" -InPlink "$WORK/plk.bed" -OutStat "$WORK/plk.ehh" -EHH 1:30000 >/dev/null 2>&1
"$NEW" -InVCF "$WORK/plk.vcf" -OutStat "$WORK/vcf.ehh" -EHH 1:30000 >/dev/null 2>&1
if cmp -s "$WORK/plk.ehh.ehh.gz" "$WORK/vcf.ehh.ehh.gz"; then
  echo "PASS  -InPlink -EHH (raw .ehh.gz bytes == -InVCF)"
else
  echo "FAIL  -InPlink -EHH"; fail=1
fi
rm -f "$WORK"/plk.ehh.* "$WORK"/vcf.ehh.*
# -OutFilterSNP: the PLINK filter listing holds the same sites as the VCF one
"$NEW" -InPlink "$WORK/plk.bed" -OutStat "$WORK/plk.flt" -OutFilterSNP >/dev/null 2>&1
"$NEW" -InVCF "$WORK/plk.vcf" -OutStat "$WORK/vcf.flt" -OutFilterSNP >/dev/null 2>&1
cut -f1,4 "$WORK/plk.flt.plink.filter.bim" > "$WORK/plk.flt.sites"
gzip -cd "$WORK/vcf.flt.vcf.filter.gz" | cut -f1,2 > "$WORK/vcf.flt.sites"
if diff "$WORK/plk.flt.sites" "$WORK/vcf.flt.sites" > "$WORK/diff.plink.flt" \
   && [ -s "$WORK/plk.flt.sites" ]; then
  echo "PASS  -InPlink -OutFilterSNP (site list == -InVCF)"
else
  echo "FAIL  -InPlink -OutFilterSNP"; head -5 "$WORK/diff.plink.flt"; fail=1
fi
# -InPlink is exclusive with the other two inputs, and a broken .bed is rejected
if "$NEW" -InPlink "$WORK/plk.bed" -InVCF "$WORK/plk.vcf" -OutStat "$WORK/plk.excl" \
     >/dev/null 2>&1; then
  echo "FAIL  -InPlink mutual exclusion"; fail=1
else
  echo "PASS  -InPlink mutual exclusion (-InVCF rejected)"
fi
printf '\x6c\x1a' > "$WORK/plkbad.bed" && cat "$WORK/plk.bed" >> "$WORK/plkbad.bed"
cp "$WORK/plk.bim" "$WORK/plkbad.bim" && cp "$WORK/plk.fam" "$WORK/plkbad.fam"
if "$NEW" -InPlink "$WORK/plkbad.bed" -OutStat "$WORK/plk.bad" >/dev/null 2>&1; then
  echo "FAIL  -InPlink individual-major .bed"; fail=1
else
  echo "PASS  -InPlink individual-major .bed rejected"
fi

# ---- real plink2 --make-bed output -------------------------------------------
# Everything above is fed by gen_sim_plink.py, which used to write a 2-byte magic
# (6c 1b) and therefore agreed with the reader's own off-by-one. No case here
# ever touched a .bed produced by the real PLINK, and that is how v0.77 shipped
# a reader that rejected every genuine PLINK dataset. These cases close that gap.
#
# What is asserted against -InVCF is what actually holds on real PLINK output:
# identical variant selection, identical NumberPairs in every distance bin, and
# r2 agreeing within a tolerance. Byte-identity does NOT hold, because a 2-bit
# .bed carries no phase: -InPlink synthesises (best, sed) per variant while
# -InVCF uses the true phase in GT. See doc/CHANGELOG.md v0.79.
PLINK2=${PLINK2:-plink2}
# Prefer real 1000 Genomes data: the simulator emits near-singleton allele counts,
# so both readers drop essentially every site and the comparison says nothing.
# With real data the two readers agree on every site and pair; see the assertions.
GT_VCF=${GT_VCF:-}
if [ -z "$GT_VCF" ] && [ -f "$CHR22_SRC" ]; then GT_VCF=$CHR22_SRC; fi
# A prepared grid VCF may also be supplied directly. This used to fall back to a
# fixed path under rebuttal/work/, which is author-side material and not part of
# this repository, so it is now an explicit override instead.
if [ -z "$GT_VCF" ] && [ -n "${GRID_VCF:-}" ] && [ -f "${GRID_VCF}" ]; then
  GT_VCF=$GRID_VCF
fi
if [ ! -x "$PLINK2" ] || [ -z "$GT_VCF" ]; then
  echo "SKIP  real plink2 --make-bed cases (plink2=$([ -x "$PLINK2" ] && echo ok || echo missing), GT_VCF=${GT_VCF:-none})"
else
    mkdir -p "$WORK/real"
    # plink2 refuses to write a .bim that holds a multiallelic variant ("p.bim
    # cannot contain multiallelic variants"), and the 1000G VCF has
    # multiallelic records. PopLDdecay2 drops those sites anyway (its Indel
    # filter), so the fixture is filtered the same way: bi-allelic SNPs only.
    # The site set the two readers are compared on is unchanged by this.
    gzip -cd "$GT_VCF" | awk 'BEGIN{FS=OFS="\t"}
        /^#/{print; next}
        { if (length($4)==1 && length($5)==1 && index($5,",")==0) print }' \
      | "$BGZIP" -c > "$WORK/real/in.vcf.gz"
    "$TABIX" -f -p vcf "$WORK/real/in.vcf.gz" 2>/dev/null
  # --vcf-half-call missing is required: a "." allele is a VCF half-call and
  # plink2 refuses to guess how to treat it.
  if "$PLINK2" --vcf-half-call missing --make-bed --vcf "$WORK/real/in.vcf.gz" \
       --out "$WORK/real/p" >/dev/null 2>&1; then
    magic=$(od -An -tx1 -N3 "$WORK/real/p.bed" | tr -d ' \n')
    nv=$(wc -l < "$WORK/real/p.bim"); ns=$(wc -l < "$WORK/real/p.fam")
    exp=$(( 3 + nv * ((ns + 3) / 4) )); act=$(stat -c%s "$WORK/real/p.bed")
    if [ "$magic" = "6c1b01" ] && [ "$exp" = "$act" ]; then
      echo "PASS  real plink2 .bed (magic $magic, $nv variants x $ns samples)"
    else
      echo "FAIL  real plink2 .bed layout (magic $magic, size $act, expected $exp)"; fail=1
    fi

    "$NEW" -InPlink "$WORK/real/p.bed" -OutStat "$WORK/real/plk" -T 1 >/dev/null 2>&1
    "$NEW" -InVCF     "$WORK/real/in.vcf.gz" -OutStat "$WORK/real/vcf" -T 1 >/dev/null 2>&1
    # NumberPairs must match bin for bin: same sites survive, same pairs counted.
    if diff <(gzip -cd "$WORK/real/plk.stat.gz" | grep -v '^#' | cut -f1,6) \
            <(gzip -cd "$WORK/real/vcf.stat.gz" | grep -v '^#' | cut -f1,6) \
            > "$WORK/real/diffpairs"; then
      echo "PASS  real plink2 .bed (NumberPairs == -InVCF, bin for bin)"
    else
      echo "FAIL  real plink2 .bed NumberPairs"; head -5 "$WORK/real/diffpairs"; fail=1
    fi
    # r2 cannot be byte-identical (a 2-bit .bed has no phase), but the pair-weighted
    # LD decay curves must agree. 1-bp bins are noise-dominated, so compare the
    # pair-weighted coarse windows, which is the measure the paper reports.
    if python3 - "$WORK/real/plk.stat.gz" "$WORK/real/vcf.stat.gz" <<'PYEOF'
import gzip, math, sys
def rd(p):
    o = []
    with gzip.open(p, 'rt') as f:
        for l in f:
            if not l.startswith('#'):
                c = l.split('\t'); o.append((int(c[0]), float(c[1]), int(c[5])))
    return o
a, b = rd(sys.argv[1]), rd(sys.argv[2])
if len(a) != len(b):
    print("  bin count differs: %d vs %d" % (len(a), len(b))); sys.exit(1)
edges = [0, 10, 50, 100, 200, 300, 600]
def agg(rows):
    s = [0.0] * (len(edges) - 1); w = [0] * (len(edges) - 1)
    for dist, r2, n in rows:
        for i in range(len(edges) - 1):
            if edges[i] <= dist < edges[i + 1]:
                s[i] += r2 * n; w[i] += n; break
    return [s[i] / w[i] for i in range(len(edges) - 1) if w[i]]
A, B = agg(a), agg(b)
if len(A) != len(B) or len(A) < 2:
    print("  coarse windows unusable"); sys.exit(1)
k = len(A); ma = sum(A) / k; mb = sum(B) / k
cov = sum((x - ma) * (y - mb) for x, y in zip(A, B))
r = cov / math.sqrt(sum((x - ma) ** 2 for x in A) * sum((y - mb) ** 2 for y in B))
print("  pair-weighted coarse-window Pearson r = %.6f" % r)
sys.exit(0 if r >= 0.999 else 1)
PYEOF
    then
      echo "PASS  real plink2 .bed (coarse-window r2 curve == -InVCF, r >= 0.999)"
    else
      echo "FAIL  real plink2 .bed r2 curve agreement"; fail=1
    fi
    # -T 8 must decode the .bed chunks exactly like -T 1
    "$NEW" -InPlink "$WORK/real/p.bed" -OutStat "$WORK/real/plk8" -T 8 >/dev/null 2>&1
    cmp_gz "real plink2 .bed (-T 8 == -T 1)" "$WORK/real/plk.stat.gz" "$WORK/real/plk8.stat.gz"
    # a mixed SNP+indel dataset (magic 6c 1b 00) must be refused, not misread
    { printf '\x6c\x1b\x00'; tail -c +4 "$WORK/real/p.bed"; } > "$WORK/real/mix.bed"
    cp "$WORK/real/p.bim" "$WORK/real/mix.bim"; cp "$WORK/real/p.fam" "$WORK/real/mix.fam"
    if "$NEW" -InPlink "$WORK/real/mix.bed" -OutStat "$WORK/real/mix" >/dev/null 2>&1; then
      echo "FAIL  -InPlink mixed SNP+indel .bed accepted"; fail=1
    else
      echo "PASS  -InPlink mixed SNP+indel .bed rejected"
    fi
    # a one-byte-short .bed must still be caught by the size check
    head -c $(( act - 1 )) "$WORK/real/p.bed" > "$WORK/real/short.bed"
    cp "$WORK/real/p.bim" "$WORK/real/short.bim"; cp "$WORK/real/p.fam" "$WORK/real/short.fam"
    if "$NEW" -InPlink "$WORK/real/short.bed" -OutStat "$WORK/real/short" >/dev/null 2>&1; then
      echo "FAIL  -InPlink truncated .bed accepted"; fail=1
    else
      echo "PASS  -InPlink truncated .bed rejected"
    fi
  else
    echo "FAIL  plink2 --make-bed failed"
    "$PLINK2" --vcf-half-call missing --make-bed --vcf "$WORK/real/in.vcf.gz" \
      --out "$WORK/real/p" 2>&1 | tr '\r' '\n' | grep -iE "error|warning" | head -5
    fail=1
  fi
fi
# Informational: -InGenotype of the same data. The genotype reader picks best/sed
# by base letter and divides MAF by 2N, the VCF/PLINK readers by observed
# alleles, so exact allele-count ties and MAF-boundary sites legitimately differ.
"$NEW" -InGenotype "$WORK/plk.geno" -OutStat "$WORK/plk.geno.stat" >/dev/null 2>&1
ndiff=$(diff <(gzip -cd "$WORK/plk.geno.stat.gz") <(gzip -cd "$WORK/plk.default.stat.gz") \
        | grep -c '^<' || true)
echo "INFO  -InGenotype vs -InPlink: $ndiff differing rows (ties/MAF boundary, expected)"
rm -f "$WORK"/plkbad.*

# ---- big-N (2504 samples, nwords=79, phased 1000G chr22 head, 10k sites) ----
# Accuracy vs old + T1/T8 consistency + SubPop. Slow to regenerate (decompresses
# the full chr22 VCF once), so results are cached under $WORK/bigN.* and only
# rebuilt when the cached marker is absent.
if [ -f "$CHR22_SRC" ]; then
  if [ ! -f "$WORK/bigN.ok" ]; then
    gzip -cd "$CHR22_SRC" \
      | awk 'BEGIN{k=0} /^#/{print; next} {k++; if(k<=10000) print; else exit}' > "$WORK/bigN.vcf"
    "$BGZIP" -f -c "$WORK/bigN.vcf" > "$WORK/bigN.vcf.gz" && rm -f "$WORK/bigN.vcf"
    "$TABIX" -p vcf "$WORK/bigN.vcf.gz"
    touch "$WORK/bigN.ok"
  fi
  "$OLD" -InVCF "$WORK/bigN.vcf.gz" -OutStat "$WORK/old.bigN.stat" 2>/dev/null
  "$NEW" -InVCF "$WORK/bigN.vcf.gz" -OutStat "$WORK/new.bigN.stat" -T 1 2>/dev/null
  cmp_gz "big-N 2504samples phased (T1 vs old)" "$WORK/old.bigN.stat.gz" "$WORK/new.bigN.stat.gz"
  "$NEW" -InVCF "$WORK/bigN.vcf.gz" -OutStat "$WORK/new.bigNt8.stat" -T 8 2>/dev/null
  cmp_gz "big-N 2504samples (T1 vs T8)" "$WORK/new.bigN.stat.gz" "$WORK/new.bigNt8.stat.gz"
  gunzip -c "$WORK/bigN.vcf.gz" | grep -m1 "^#CHROM" | cut -f10- | tr '\t' '\n' | head -500 | tr '\n' ' ' > "$WORK/bigN.g1.list"
  "$OLD" -InVCF "$WORK/bigN.vcf.gz" -OutStat "$WORK/old.bigNsub.stat" -SubPop "$WORK/bigN.g1.list" 2>/dev/null
  "$NEW" -InVCF "$WORK/bigN.vcf.gz" -OutStat "$WORK/new.bigNsub.stat" -SubPop "$WORK/bigN.g1.list" -T 1 2>/dev/null
  cmp_gz "big-N 2504samples -SubPop 500 (T1 vs old)" "$WORK/old.bigNsub.stat.gz" "$WORK/new.bigNsub.stat.gz"
  # -InPlink on the same big-N data: the .bed/.bim/.fam conversion keeps only the
  # records PLINK can express (9723 of the 10000), and writes the matching VCF.
  # -InPlink must equal -InVCF on that VCF, which in turn equals 3.45.
  if [ ! -f "$WORK/bigNplk.bed" ]; then
    python3 RunTest/06.Verify/gen_sim_plink.py -o "$WORK/bigNplk" -i "$WORK/bigN.vcf.gz" >/dev/null
  fi
  cut -f2 "$WORK/bigNplk.fam" | head -500 | tr '\n' ' ' > "$WORK/bigNplk.g1.list"
  "$OLD" -InVCF "$WORK/bigNplk.vcf" -OutStat "$WORK/old.bigNplk.stat" 2>/dev/null
  "$NEW" -InVCF "$WORK/bigNplk.vcf" -OutStat "$WORK/new.bigNplkvc.stat" -T 1 2>/dev/null
  "$NEW" -InPlink "$WORK/bigNplk.bed" -OutStat "$WORK/new.bigNplk.stat" -T 1 2>/dev/null
  if diff <(gzip -cd "$WORK/new.bigNplk.stat.gz") <(gzip -cd "$WORK/new.bigNplkvc.stat.gz") > "$WORK/diff.bigNplk" \
     && diff <(gzip -cd "$WORK/new.bigNplk.stat.gz") <(gzip -cd "$WORK/old.bigNplk.stat.gz") >/dev/null; then
    echo "PASS  big-N -InPlink 2504samples (== -InVCF == old)"
  else
    echo "FAIL  big-N -InPlink"; head -20 "$WORK/diff.bigNplk"; fail=1
  fi
  "$NEW" -InPlink "$WORK/bigNplk.bed" -OutStat "$WORK/new.bigNplkt8.stat" -T 8 2>/dev/null
  cmp_gz "big-N -InPlink 2504samples (T1 vs T8)" "$WORK/new.bigNplk.stat.gz" "$WORK/new.bigNplkt8.stat.gz"
  "$NEW" -InPlink "$WORK/bigNplk.bed" -OutStat "$WORK/new.bigNplksub.stat" -SubPop "$WORK/bigNplk.g1.list" -T 1 2>/dev/null
  "$NEW" -InVCF "$WORK/bigNplk.vcf" -OutStat "$WORK/new.bigNplkvcsub.stat" -SubPop "$WORK/bigNplk.g1.list" -T 1 2>/dev/null
  cmp_gz "big-N -InPlink 2504samples -SubPop 500 (== -InVCF)" "$WORK/new.bigNplksub.stat.gz" "$WORK/new.bigNplkvcsub.stat.gz"
else
  echo "SKIP  big-N 2504samples (chr22 data not present)"
fi

exit $fail