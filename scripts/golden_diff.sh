#!/bin/bash
# Golden-diff: PopLDdecay2 output must be byte-identical (after gunzip) to the
# original PopLDdecay for the same input / parameters.
set -uo pipefail
cd "$(dirname "$0")/.."

OLD=PopLDdecay-3.45/bin/PopLDdecay
NEW=bin/PopLDdecay2
BGZIP=/home/heweiming/01.Software/samtools-1.23/bin/bgzip
TABIX=/home/heweiming/01.Software/samtools-1.23/bin/tabix
WORK=Test//ldgolden
mkdir -p "$WORK"

fail=0

# ---- inputs (plain / gzip / bgzip / phased) ----
python3 scripts/gen_sim_vcf.py -o "$WORK/sim.vcf" --nsamples 60 --sites 2000 --chroms 4 >/dev/null
python3 scripts/gen_sim_vcf.py -o "$WORK/sim.vcf.gz" --gzip --nsamples 60 --sites 2000 --chroms 4 >/dev/null
"$BGZIP" -f -c "$WORK/sim.vcf" > "$WORK/sim.bgzip.vcf.gz"
python3 scripts/gen_sim_vcf.py -o "$WORK/sim_phased.vcf.gz" --gzip --phased --nsamples 60 --sites 2000 --chroms 4 >/dev/null

for name in sim.vcf sim.vcf.gz sim.bgzip.vcf.gz sim_phased.vcf.gz; do
  in="$WORK/$name"
  "$OLD" -InVCF "$in" -OutStat "$WORK/old.$name.stat" 2>/dev/null
  "$NEW" -InVCF "$in" -OutStat "$WORK/new.$name.stat" 2>/dev/null
  if diff <(gzip -cd "$WORK/old.$name.stat.gz") <(gzip -cd "$WORK/new.$name.stat.gz") > "$WORK/diff.$name"; then
    echo "PASS  $name"
  else
    echo "FAIL  $name"
    head -20 "$WORK/diff.$name"
    fail=1
  fi
done

# ---- parameter variation ----
"$OLD" -InVCF "$WORK/sim.vcf.gz" -OutStat "$WORK/old.param.stat" -MAF 0.1 -Het 0.5 -Miss 0.1 -MaxDist 100 2>/dev/null
"$NEW" -InVCF "$WORK/sim.vcf.gz" -OutStat "$WORK/new.param.stat" -MAF 0.1 -Het 0.5 -Miss 0.1 -MaxDist 100 2>/dev/null
if diff <(gzip -cd "$WORK/old.param.stat.gz") <(gzip -cd "$WORK/new.param.stat.gz") > "$WORK/diff.param"; then
  echo "PASS  param(-MAF 0.1 -Het 0.5 -Miss 0.1 -MaxDist 100)"
else
  echo "FAIL  param"; head -20 "$WORK/diff.param"; fail=1
fi

# ---- default MaxDist (300kb) is covered by sim runs; also test -OutStat already ending in .stat.gz ----
"$OLD" -InVCF "$WORK/sim.vcf.gz" -OutStat "$WORK/old.prefix.stat.gz" 2>/dev/null
"$NEW" -InVCF "$WORK/sim.vcf.gz" -OutStat "$WORK/new.prefix.stat.gz" 2>/dev/null
if diff <(gzip -cd "$WORK/old.prefix.stat.gz") <(gzip -cd "$WORK/new.prefix.stat.gz") > "$WORK/diff.prefix"; then
  echo "PASS  prefix(-OutStat *.stat.gz)"
else
  echo "FAIL  prefix"; head -20 "$WORK/diff.prefix"; fail=1
fi

# ---- multithreaded self-consistency: -T 8 must equal -T 1 (formatted) ----
"$NEW" -InVCF "$WORK/sim.vcf.gz" -OutStat "$WORK/t1.stat" -T 1 2>/dev/null
"$NEW" -InVCF "$WORK/sim.vcf.gz" -OutStat "$WORK/t8.stat" -T 8 2>/dev/null
if diff <(gzip -cd "$WORK/t1.stat.gz") <(gzip -cd "$WORK/t8.stat.gz") > "$WORK/diff.threads"; then
  echo "PASS  threads(-T 8 vs -T 1)"
else
  echo "FAIL  threads"; head -20 "$WORK/diff.threads"; fail=1
fi

# ---- indexed fast path: bgzip + tabix index, -T 4 must match old (fast path) ----
"$BGZIP" -f -c "$WORK/sim.vcf" > "$WORK/sim.idx.vcf.gz"
"$TABIX" -p vcf "$WORK/sim.idx.vcf.gz"
"$OLD" -InVCF "$WORK/sim.idx.vcf.gz" -OutStat "$WORK/old.simidx.stat" 2>/dev/null
"$NEW" -InVCF "$WORK/sim.idx.vcf.gz" -OutStat "$WORK/new.simidx.stat" -T 4 2>/dev/null
if diff <(gzip -cd "$WORK/old.simidx.stat.gz") <(gzip -cd "$WORK/new.simidx.stat.gz") > "$WORK/diff.simidx"; then
  echo "PASS  fastpath indexed(bgzip+tabix, -T 4 vs old)"
else
  echo "FAIL  fastpath indexed"; head -20 "$WORK/diff.simidx"; fail=1
fi
# sharded fastpath: threads(8) > contigs(4) -> position shards; must still match old
"$OLD" -InVCF "$WORK/sim.idx.vcf.gz" -OutStat "$WORK/old.shard.stat" 2>/dev/null
"$NEW" -InVCF "$WORK/sim.idx.vcf.gz" -OutStat "$WORK/new.shard.stat" -T 8 2>/dev/null
if diff <(gzip -cd "$WORK/old.shard.stat.gz") <(gzip -cd "$WORK/new.shard.stat.gz") > "$WORK/diff.shard"; then
  echo "PASS  fastpath sharded(threads>contigs, -T 8 vs old)"
else
  echo "FAIL  fastpath sharded"; head -20 "$WORK/diff.shard"; fail=1
fi
# phased indexed counterpart (for the fastpath .LD.gz checks below)
gzip -cd "$WORK/sim_phased.vcf.gz" | "$BGZIP" -f -c > "$WORK/sim_phased.idx.vcf.gz"
"$TABIX" -p vcf "$WORK/sim_phased.idx.vcf.gz"
# sharded fastpath .LD.gz byte-identical to the fallback (OutType 7, phased)
"$NEW" -InVCF "$WORK/sim_phased.idx.vcf.gz" -OutStat "$WORK/new.shardld.t1.stat" -OutType 7 -T 1 2>/dev/null
"$NEW" -InVCF "$WORK/sim_phased.idx.vcf.gz" -OutStat "$WORK/new.shardld.t8.stat" -OutType 7 -T 8 2>/dev/null
if diff <(gzip -cd "$WORK/new.shardld.t1.LD.gz") <(gzip -cd "$WORK/new.shardld.t8.LD.gz") > "$WORK/diff.shardld"; then
  echo "PASS  fastpath sharded (OutType 7 .LD.gz T8 vs T1)"
else
  echo "FAIL  fastpath sharded .LD.gz"; head -20 "$WORK/diff.shardld"; fail=1
fi
# single-chromosome indexed fastpath (0.6.8): contigs>=1 now admitted; threads
# (8) > 1 chr -> position shards on the single chromosome; must match old
python3 scripts/gen_sim_vcf.py -o "$WORK/sim1.vcf" --nsamples 60 --sites 1500 --chroms 1 >/dev/null
"$BGZIP" -f -c "$WORK/sim1.vcf" > "$WORK/sim1.idx.vcf.gz"
"$TABIX" -p vcf "$WORK/sim1.idx.vcf.gz"
"$OLD" -InVCF "$WORK/sim1.idx.vcf.gz" -OutStat "$WORK/old.sim1.stat" 2>/dev/null
"$NEW" -InVCF "$WORK/sim1.idx.vcf.gz" -OutStat "$WORK/new.sim1.stat" -T 8 2>/dev/null
if diff <(gzip -cd "$WORK/old.sim1.stat.gz") <(gzip -cd "$WORK/new.sim1.stat.gz") > "$WORK/diff.sim1"; then
  echo "PASS  fastpath single-chr (bgzip+tabix, -T 8 vs old)"
else
  echo "FAIL  fastpath single-chr"; head -20 "$WORK/diff.sim1"; fail=1
fi

# ---- OutType 0-8: .stat.gz (+ .LD.gz for 3/6/7/8) must match old, unphased & phased ----
for ot in 0 1 2 3 4 5 6 7 8; do
  for phase in unph phased; do
    if [ "$phase" = phased ]; then in="$WORK/sim_phased.vcf.gz"; else in="$WORK/sim.vcf.gz"; fi
    "$OLD" -InVCF "$in" -OutStat "$WORK/old.ot$ot.$phase.stat" -OutType "$ot" 2>/dev/null
    "$NEW" -InVCF "$in" -OutStat "$WORK/new.ot$ot.$phase.stat" -OutType "$ot" -T 1 2>/dev/null
    if diff <(gzip -cd "$WORK/old.ot$ot.$phase.stat.gz") <(gzip -cd "$WORK/new.ot$ot.$phase.stat.gz") > "$WORK/diff.ot$ot.$phase"; then
      echo "PASS  OutType $ot $phase (.stat.gz)"
    else
      echo "FAIL  OutType $ot $phase (.stat.gz)"; head -20 "$WORK/diff.ot$ot.$phase"; fail=1
    fi
    case "$ot" in
      3|6|7|8)
        if diff <(gzip -cd "$WORK/old.ot$ot.$phase.LD.gz") <(gzip -cd "$WORK/new.ot$ot.$phase.LD.gz") > "$WORK/diff.ld.ot$ot.$phase"; then
          echo "PASS  OutType $ot $phase (.LD.gz)"
        else
          echo "FAIL  OutType $ot $phase (.LD.gz)"; head -20 "$WORK/diff.ld.ot$ot.$phase"; fail=1
        fi
        # indexed fast path (-T 4) must produce the identical .LD.gz
        if [ "$phase" = phased ]; then idxin="$WORK/sim_phased.idx.vcf.gz"; else idxin="$WORK/sim.idx.vcf.gz"; fi
        "$NEW" -InVCF "$idxin" -OutStat "$WORK/new.ot$ot.fast.stat" -OutType "$ot" -T 4 2>/dev/null
        if diff <(gzip -cd "$WORK/new.ot$ot.$phase.LD.gz") <(gzip -cd "$WORK/new.ot$ot.fast.LD.gz") > "$WORK/diff.ld.fast.ot$ot"; then
          echo "PASS  OutType $ot $phase (fastpath .LD.gz)"
        else
          echo "FAIL  OutType $ot $phase (fastpath .LD.gz)"; head -20 "$WORK/diff.ld.fast.ot$ot"; fail=1
        fi
        ;;
    esac
  done
done

# ---- -L region restrict (new feature, no old-binary target: self-consistency) ----
# (a) fallback vs fastpath on the indexed input must be formatted-identical
"$NEW" -InVCF "$WORK/sim.idx.vcf.gz" -OutStat "$WORK/l1.t1.stat" -L chr1 -T 1 2>/dev/null
"$NEW" -InVCF "$WORK/sim.idx.vcf.gz" -OutStat "$WORK/l1.t4.stat" -L chr1 -T 4 2>/dev/null
if diff <(gzip -cd "$WORK/l1.t1.stat.gz") <(gzip -cd "$WORK/l1.t4.stat.gz") > "$WORK/diff.l1"; then
  echo "PASS  -L chr1 (fallback T1 vs fastpath T4)"
else
  echo "FAIL  -L chr1"; head -20 "$WORK/diff.l1"; fail=1
fi
# (b) -L chr1 on the full VCF == running on a chr1-only sub-VCF (byte-identical)
awk -F'\t' '$1=="#CHROM" || $1=="chr1"' "$WORK/sim.vcf" > "$WORK/sim.chr1.vcf"
"$NEW" -InVCF "$WORK/sim.vcf" -OutStat "$WORK/eq.full.stat" -L chr1 -T 1 2>/dev/null
"$NEW" -InVCF "$WORK/sim.chr1.vcf" -OutStat "$WORK/eq.sub.stat" -T 1 2>/dev/null
if diff <(gzip -cd "$WORK/eq.full.stat.gz") <(gzip -cd "$WORK/eq.sub.stat.gz") > "$WORK/diff.eq"; then
  echo "PASS  -L chr1 == chr1-only sub-VCF"
else
  echo "FAIL  -L chr1 == sub-VCF"; head -20 "$WORK/diff.eq"; fail=1
fi
# (c) BED file (0-based half-open) fallback vs fastpath
printf 'chr1\t999\t1000000\nchr2\t1\t5000000\n' > "$WORK/regions.bed"
"$NEW" -InVCF "$WORK/sim.idx.vcf.gz" -OutStat "$WORK/bed.t1.stat" -L "$WORK/regions.bed" -T 1 2>/dev/null
"$NEW" -InVCF "$WORK/sim.idx.vcf.gz" -OutStat "$WORK/bed.t4.stat" -L "$WORK/regions.bed" -T 4 2>/dev/null
if diff <(gzip -cd "$WORK/bed.t1.stat.gz") <(gzip -cd "$WORK/bed.t4.stat.gz") > "$WORK/diff.bed"; then
  echo "PASS  -L BED (fallback T1 vs fastpath T4)"
else
  echo "FAIL  -L BED"; head -20 "$WORK/diff.bed"; fail=1
fi
# (d) empty region (chr not present) -> header-only output, both paths equal
"$NEW" -InVCF "$WORK/sim.idx.vcf.gz" -OutStat "$WORK/empty.t1.stat" -L chrY -T 1 2>/dev/null
"$NEW" -InVCF "$WORK/sim.idx.vcf.gz" -OutStat "$WORK/empty.t4.stat" -L chrY -T 4 2>/dev/null
if diff <(gzip -cd "$WORK/empty.t1.stat.gz") <(gzip -cd "$WORK/empty.t4.stat.gz") > "$WORK/diff.empty"; then
  echo "PASS  -L empty region (header only)"
else
  echo "FAIL  -L empty region"; head -20 "$WORK/diff.empty"; fail=1
fi

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
if diff <(gzip -cd "$WORK/dupreg.t1.stat.gz") <(gzip -cd "$WORK/dupreg.t4.stat.gz") > "$WORK/diff.dupreg"; then
  echo "PASS  duplicate -L chr1,chr1 (T1 vs T4)"
else
  echo "FAIL  duplicate -L chr1,chr1"; head -20 "$WORK/diff.dupreg"; fail=1
fi

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
if diff <(gzip -cd "$WORK/haploid_mid.old.stat.gz") <(gzip -cd "$WORK/haploid_mid.new.stat.gz") > "$WORK/diff.haploid_mid"; then
  echo "PASS  single-char GT (mid columns, vs old)"
else
  echo "FAIL  single-char GT mid columns"; head -20 "$WORK/diff.haploid_mid"; fail=1
fi

# ---- -SubPop subset (whitespace-separated list, 30 of 60 samples) ----
{ for i in $(seq 1 2 59); do echo -n "S$i "; done; echo; echo "# comment"; echo "S99"; } > "$WORK/pop1.list"
for phase in unphased phased; do
  if [ "$phase" = phased ]; then in="$WORK/sim_phased.vcf.gz"; else in="$WORK/sim.vcf.gz"; fi
  "$OLD" -InVCF "$in" -OutStat "$WORK/old.sub.$phase.stat" -SubPop "$WORK/pop1.list" 2>/dev/null
  "$NEW" -InVCF "$in" -OutStat "$WORK/new.sub.$phase.stat" -SubPop "$WORK/pop1.list" -T 1 2>/dev/null
  if diff <(gzip -cd "$WORK/old.sub.$phase.stat.gz") <(gzip -cd "$WORK/new.sub.$phase.stat.gz") > "$WORK/diff.sub.$phase"; then
    echo "PASS  -SubPop $phase (30/60 vs old)"
  else
    echo "FAIL  -SubPop $phase"; head -20 "$WORK/diff.sub.$phase"; fail=1
  fi
done
# SubPop on the indexed input: fastpath (T4) must equal fallback (T1)
"$NEW" -InVCF "$WORK/sim.idx.vcf.gz" -OutStat "$WORK/sub.t1.stat" -SubPop "$WORK/pop1.list" -T 1 2>/dev/null
"$NEW" -InVCF "$WORK/sim.idx.vcf.gz" -OutStat "$WORK/sub.t4.stat" -SubPop "$WORK/pop1.list" -T 4 2>/dev/null
if diff <(gzip -cd "$WORK/sub.t1.stat.gz") <(gzip -cd "$WORK/sub.t4.stat.gz") > "$WORK/diff.sub.t"; then
  echo "PASS  -SubPop indexed (fallback T1 vs fastpath T4)"
else
  echo "FAIL  -SubPop indexed"; head -20 "$WORK/diff.sub.t"; fail=1
fi
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
CHR22_SRC=Test/1000genomes_chr22/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
if [ -f "$CHR22_SRC" ]; then
  if [ ! -f "$WORK/chr22head.vcf.gz" ]; then
    gzip -cd "$CHR22_SRC" | awk 'BEGIN{k=0} /^#/{print; next} {k++; if(k<=2000) print; else exit}' > "$WORK/chr22head.vcf"
    "$BGZIP" -f -c "$WORK/chr22head.vcf" > "$WORK/chr22head.vcf.gz" && rm -f "$WORK/chr22head.vcf"
    "$TABIX" -p vcf "$WORK/chr22head.vcf.gz"
  fi
  for T in 1 4 8; do
    "$OLD" -InVCF "$WORK/chr22head.vcf.gz" -OutStat "$WORK/old.c22.t$T.stat" 2>/dev/null
    "$NEW" -InVCF "$WORK/chr22head.vcf.gz" -OutStat "$WORK/new.c22.t$T.stat" -T "$T" 2>/dev/null
    if diff <(gzip -cd "$WORK/old.c22.t$T.stat.gz") <(gzip -cd "$WORK/new.c22.t$T.stat.gz") > "$WORK/diff.c22.t$T"; then
      echo "PASS  real chr22 head -T $T"
    else
      echo "FAIL  real chr22 head -T $T"; head -20 "$WORK/diff.c22.t$T"; fail=1
    fi
  done
  # SubPop: first 500 of the 2504 samples (matches old; fastpath == fallback)
  gunzip -c "$CHR22_SRC" | grep -m1 "^#CHROM" | cut -f10- | tr '\t' '\n' \
    | head -500 | tr '\n' ' ' > "$WORK/c22.g1.list"
  "$OLD" -InVCF "$WORK/chr22head.vcf.gz" -OutStat "$WORK/old.c22.sub.stat" -SubPop "$WORK/c22.g1.list" 2>/dev/null
  "$NEW" -InVCF "$WORK/chr22head.vcf.gz" -OutStat "$WORK/new.c22.sub.stat" -SubPop "$WORK/c22.g1.list" -T 1 2>/dev/null
  if diff <(gzip -cd "$WORK/old.c22.sub.stat.gz") <(gzip -cd "$WORK/new.c22.sub.stat.gz") > "$WORK/diff.c22.sub"; then
    echo "PASS  real chr22 head -SubPop (500 samples vs old)"
  else
    echo "FAIL  real chr22 head -SubPop"; head -20 "$WORK/diff.c22.sub"; fail=1
  fi
  "$NEW" -InVCF "$WORK/chr22head.vcf.gz" -OutStat "$WORK/new.c22.subt.stat" -SubPop "$WORK/c22.g1.list" -T 8 2>/dev/null
  if diff <(gzip -cd "$WORK/new.c22.sub.stat.gz") <(gzip -cd "$WORK/new.c22.subt.stat.gz") > "$WORK/diff.c22.subt"; then
    echo "PASS  real chr22 head -SubPop indexed (fallback T1 vs fastpath T8)"
  else
    echo "FAIL  real chr22 head -SubPop indexed"; head -20 "$WORK/diff.c22.subt"; fail=1
  fi
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

# ---- -InGenotype (M6): stat content must match old; filter.gz raw bytes ----
geno_ok() {  # name infile extra...
  local name="$1" in="$2"; shift 2
  "$OLD" -InGenotype "$in" -OutStat "$WORK/old.$name.stat" "$@" >/dev/null 2>&1
  "$NEW" -InGenotype "$in" -OutStat "$WORK/new.$name.stat" "$@" >/dev/null 2>&1
  if diff <(gzip -cd "$WORK/old.$name.stat.gz") <(gzip -cd "$WORK/new.$name.stat.gz") > "$WORK/diff.$name"; then
    echo "PASS  -InGenotype $name (.stat.gz == old)"
  else
    echo "FAIL  -InGenotype $name"; head -20 "$WORK/diff.$name"; fail=1
  fi
}
python3 scripts/gen_sim_genotype.py -i "$WORK/sim.vcf.gz" -o "$WORK/sim.geno" --with-header >/dev/null
geno_ok geno "$WORK/sim.geno"
# gz input (plain-gzip genotype input)
gzip -c "$WORK/sim.geno" > "$WORK/sim.geno.gz"
geno_ok geno.gz "$WORK/sim.geno.gz"
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
if diff <(gzip -cd "$WORK/old.genomix.stat.gz") <(gzip -cd "$WORK/new.genomix.stat.gz") > "$WORK/diff.genomix"; then
  echo "PASS  -InGenotype mixed with -InVCF (.stat.gz == old)"
else
  echo "FAIL  -InGenotype mixed with -InVCF"; head -20 "$WORK/diff.genomix"; fail=1
fi
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
  if diff <(gzip -cd "$WORK/old.bigN.stat.gz") <(gzip -cd "$WORK/new.bigN.stat.gz") > "$WORK/diff.bigN"; then
    echo "PASS  big-N 2504samples phased (T1 vs old)"
  else
    echo "FAIL  big-N T1"; head -20 "$WORK/diff.bigN"; fail=1
  fi
  "$NEW" -InVCF "$WORK/bigN.vcf.gz" -OutStat "$WORK/new.bigNt8.stat" -T 8 2>/dev/null
  if diff <(gzip -cd "$WORK/new.bigN.stat.gz") <(gzip -cd "$WORK/new.bigNt8.stat.gz") > "$WORK/diff.bigN.t"; then
    echo "PASS  big-N 2504samples (T1 vs T8)"
  else
    echo "FAIL  big-N T8"; head -20 "$WORK/diff.bigN.t"; fail=1
  fi
  gunzip -c "$WORK/bigN.vcf.gz" | grep -m1 "^#CHROM" | cut -f10- | tr '\t' '\n' | head -500 | tr '\n' ' ' > "$WORK/bigN.g1.list"
  "$OLD" -InVCF "$WORK/bigN.vcf.gz" -OutStat "$WORK/old.bigNsub.stat" -SubPop "$WORK/bigN.g1.list" 2>/dev/null
  "$NEW" -InVCF "$WORK/bigN.vcf.gz" -OutStat "$WORK/new.bigNsub.stat" -SubPop "$WORK/bigN.g1.list" -T 1 2>/dev/null
  if diff <(gzip -cd "$WORK/old.bigNsub.stat.gz") <(gzip -cd "$WORK/new.bigNsub.stat.gz") > "$WORK/diff.bigNsub"; then
    echo "PASS  big-N 2504samples -SubPop 500 (T1 vs old)"
  else
    echo "FAIL  big-N -SubPop"; head -20 "$WORK/diff.bigNsub"; fail=1
  fi
else
  echo "SKIP  big-N 2504samples (chr22 data not present)"
fi

exit $fail