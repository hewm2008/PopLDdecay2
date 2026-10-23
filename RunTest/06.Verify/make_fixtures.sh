#!/bin/bash
# Regenerate expected.tsv: freeze what the current binary produces for every
# case, together with the exact command line. step1_checkOutput.sh then rebuilds the same
# inputs and checks its output against these digests.
#
#   bash make_fixtures.sh   # regenerate (needs this tree)
#   bash step1_checkOutput.sh             # what CI runs
#
# Format (tab separated): case_id, kind (stat|raw), argv, output_suffix, sha256
#   stat = digest of `zcat <out>`; raw = digest of the file itself
set -uo pipefail
# This file lives in <repo>/RunTest/06.Verify/ -> the repository root is ../..
# Everything it needs is derived here or taken from PATH; it deliberately does
# NOT source the rebuttal environment, which hardcodes absolute paths to a
# non-redistributable PopLDdecay 3.45 build and to third-party tools.
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE/../.."
REPO=$(pwd)
NEW=$REPO/bin/PopLDdecay2
WORK=${WORK:-$(mktemp -d)}
PY=${PYTHON:-python3}
BGZIP=${BGZIP:-bgzip}
TABIX=${TABIX:-tabix}
for _t in "$PY" "$BGZIP" "$TABIX"; do
  command -v "$_t" >/dev/null 2>&1 || {
    echo "FATAL: '$_t' not found in PATH."
    echo "       Install python3 and samtools (which provides bgzip and tabix),"
    echo "       or point the matching variable at the binary and re-run."
    exit 2
  }
done
# Resolve to absolute paths: gen_inputs.sh builds several fixtures only when
# [ -f "$BGZIP" ], which is false for a bare command name.
BGZIP=$(command -v "$BGZIP")
TABIX=$(command -v "$TABIX")
export REPO NEW WORK PY BGZIP TABIX
# Input generation lives in one place and is shared with the CI check, so the
# inputs a digest was computed from are exactly the inputs CI reconstructs.
. "$HERE/gen_inputs.sh"

FIXDIR=$HERE
WORKD=$WORK/fixtures
mkdir -p "$FIXDIR" "$WORKD"
[ -x "$NEW" ] || { echo "build bin/PopLDdecay2 first" >&2; exit 1; }

sha_stat() { zcat "$1" 2>/dev/null | sha256sum | cut -d' ' -f1; }
sha_raw()  { sha256sum "$1" 2>/dev/null | cut -d' ' -f1; }

emit() {   # id kind argv produced_file suffix
  local id=$1 kind=$2 argv=$3 file=$4 suffix=$5
  if [ ! -f "$file" ]; then echo "  SKIP $id (missing $file)" >&2; return 0; fi
  local h
  [ "$kind" = stat ] && h=$(sha_stat "$file") || h=$(sha_raw "$file")
  printf '%s\t%s\t%s\t%s\t%s\n' "$id" "$kind" "$argv" "$suffix" "$h" >> "$OUT"
  printf '  %-30s %-20s %s\n' "$id" "$suffix" "${h:0:16}"
}

# ---------------------------------------------------------------- inputs
make_ci_inputs "$WORKD"
V=$WORKD/sim.vcf.gz
VP=$WORKD/sim_phased.vcf.gz
VG=$WORKD/sim.bgz.vcf.gz
GT=$WORKD/sim.geno.txt
SUB=$WORKD/sub30.txt
GRP=$WORKD/groups2.txt

[ -f "$WORKD/sim.vcf" ] || python3 "$REPO/RunTest/06.Verify/gen_sim_vcf.py" -o "$WORKD/sim.vcf" \
                         --nsamples 60 --sites 2000 --chroms 4 >/dev/null 2>&1
[ -f "$V" ]  || gzip -c "$WORKD/sim.vcf" > "$V"
[ -f "$VP" ] || python3 "$REPO/RunTest/06.Verify/gen_sim_vcf.py" -o "$WORKD/sim_phased.vcf.gz" --gzip \
                         --phased --nsamples 60 --sites 2000 --chroms 4 >/dev/null 2>&1
[ -f "$VG" ] || "$BGZIP" -f -c "$WORKD/sim.vcf" > "$VG"
[ -f "$VG.tbi" ] || "$TABIX" -f -p vcf "$VG"
[ -f "$GT" ] || python3 "$REPO/RunTest/06.Verify/gen_sim_genotype.py" -i "$WORKD/sim.vcf" \
                        -o "$GT" --nsamples 60 >/dev/null 2>&1
# The generator keeps every sample on the single #CHROM line: split them out.
samps() { zcat "$V" 2>/dev/null | grep -m1 '^#CHROM' | cut -f10- | tr '\t' '\n'; }
[ -f "$SUB" ] || samps | head -30 > "$SUB"
[ -f "$GRP" ] || samps | head -40 | awk '{print $1"\t"((NR<=20)?"G1":"G2")}' > "$GRP"

# ---------------------------------------------------------------- expected
# Case families mirror golden_diff.sh one for one, so the matrix the
# manuscript describes has exactly one counterpart per row here:
#   golden_diff.sh  compares PopLDdecay2 against PopLDdecay 3.45 (needs the
#                   non-redistributable 3.45 binary, run by the authors)
#   this file       freezes PopLDdecay2's own output for the same command
#                   lines, so CI can run the identical matrix with no 3.45
# The two share the argv of every case; only the reference differs.
OUT=$FIXDIR/expected.tsv
{
printf '# PopLDdecay2 expected-output fixtures -- CI regression gate.\n'
printf '#\n'
printf '# Regenerate: bash RunTest/06.Verify/make_fixtures.sh\n'
printf '# Verify   : bash RunTest/06.Verify/step1_checkOutput.sh   (CI; no 3.45 needed)\n'
printf '#\n'
printf '# PopLDdecay 3.45 is not redistributable, so CI cannot re-run it. These\n'
printf '# digests come from the current implementation, whose byte-for-byte equality\n'
printf '# with PopLDdecay 3.45 over the same case matrix is verified separately by\n'
printf '# golden_diff.sh. Regenerate the fixtures whenever an output format\n'
printf '# changes *intentionally*, and record why in the commit message.\n'
printf '#\n'
printf '# case_id\tkind\targv\toutput_suffix\tsha256\n'
} > "$OUT"

echo "=== generating fixtures ==="

# ---- 1. input containers ----
for src in sim.vcf sim.vcf.gz sim.bgz.vcf.gz sim_phased.vcf.gz; do
  [ -f "$WORKD/$src" ] || continue
  o=$WORKD/out_in_$(echo "$src" | tr '.' '_')
  "$NEW" -InVCF "$WORKD/$src" -OutStat "$o" -MaxDist 100 >/dev/null 2>&1
  emit "input/$src" stat "-InVCF @WORK@/$src -OutStat @OUT@ -MaxDist 100" "$o.stat.gz" stat.gz
done

# ---- 2. filter parameters ----
param_case() {   # id  key=value ...
  local id=$1; shift
  local o=$WORKD/out_parm_$(echo "$id" | tr '/:.' '___')
  local -a a=(); local kv f=""
  for kv in "$@"; do a+=("-${kv%%=*}" "${kv#*=}"); f="$f -${kv%%=*} ${kv#*=}"; done
  "$NEW" -InVCF "$V" -OutStat "$o" "${a[@]}" >/dev/null 2>&1
  emit "param/$id" stat "-InVCF @WORK@/sim.vcf.gz -OutStat @OUT@$f" "$o.stat.gz" stat.gz
}
param_case MAF_0.1     MAF=0.1
param_case Het_0.5     Het=0.5
param_case Miss_0.1    Miss=0.1
param_case MaxDist_50  MaxDist=50
param_case combo       MAF=0.05 Het=0.8 Miss=0.2 MaxDist=200

# ---- 3. OutType 0-8 ----
for ot in 0 1 2 3 4 5 6 7 8; do
  o=$WORKD/out.ot$ot
  "$NEW" -InVCF "$V" -OutStat "$o" -OutType $ot -MaxDist 100 >/dev/null 2>&1
  emit "outtype/$ot/stat" stat "-InVCF @WORK@/sim.vcf.gz -OutStat @OUT@ -OutType $ot -MaxDist 100" "$o.stat.gz" stat.gz
done

# ---- 4. OutType 3/6/7/8 pairwise table ----
for ot in 3 6 7 8; do
  o=$WORKD/out.ot${ot}p
  "$NEW" -InVCF "$V" -OutStat "$o" -OutType $ot -MaxDist 50 >/dev/null 2>&1
  emit "outtype/$ot/LD" raw "-InVCF @WORK@/sim.vcf.gz -OutStat @OUT@ -OutType $ot -MaxDist 50" "$o.LD.gz" LD.gz
done

# ---- 5. phased input ----
for ot in 1 2 3 4 8; do
  o=$WORKD/out.ph$ot
  "$NEW" -InVCF "$VP" -OutStat "$o" -OutType $ot -MaxDist 100 >/dev/null 2>&1
  emit "phased/outtype/$ot" stat "-InVCF @WORK@/sim_phased.vcf.gz -OutStat @OUT@ -OutType $ot -MaxDist 100" "$o.stat.gz" stat.gz
done

# ---- 6. thread consistency: these four MUST digest identically ----
for t in 2 4 8; do
  o=$WORKD/out.t$t
  "$NEW" -InVCF "$VG" -OutStat "$o" -T $t -MaxDist 100 >/dev/null 2>&1
  emit "threads/$t" stat "-InVCF @WORK@/sim.bgz.vcf.gz -OutStat @OUT@ -T $t -MaxDist 100" "$o.stat.gz" stat.gz
done

# ---- 7. region restriction ----
"$NEW" -InVCF "$VG" -OutStat "$WORKD/out.L1" -L chr1        -T 4 -MaxDist 100 >/dev/null 2>&1
emit "region/chr1"       stat "-InVCF @WORK@/sim.bgz.vcf.gz -OutStat @OUT@ -L chr1 -T 4 -MaxDist 100"        "$WORKD/out.L1.stat.gz" stat.gz
"$NEW" -InVCF "$VG" -OutStat "$WORKD/out.L2" -L chr1:1-400000 -T 4 -MaxDist 100 >/dev/null 2>&1
emit "region/chr1_range" stat "-InVCF @WORK@/sim.bgz.vcf.gz -OutStat @OUT@ -L chr1:1-400000 -T 4 -MaxDist 100" "$WORKD/out.L2.stat.gz" stat.gz

# ---- 8. SubPop (single column keeps the legacy name; two column labels groups) ----
o=$WORKD/out_sub
"$NEW" -InVCF "$V" -OutStat "$o" -SubPop "$SUB" -MaxDist 100 >/dev/null 2>&1
emit "subpop/30of60"    stat "-InVCF @WORK@/sim.vcf.gz -OutStat @OUT@ -SubPop @WORK@/sub30.txt -MaxDist 100"  "$o.stat.gz"    stat.gz
o2=$WORKD/out_sub2
"$NEW" -InVCF "$V" -OutStat "$o2" -SubPop "$GRP" -MaxDist 100 >/dev/null 2>&1
emit "subpop/twocol/G1" stat "-InVCF @WORK@/sim.vcf.gz -OutStat @OUT@ -SubPop @WORK@/groups2.txt -MaxDist 100" "$o2.G1.stat.gz" G1.stat.gz
emit "subpop/twocol/G2" stat "-InVCF @WORK@/sim.vcf.gz -OutStat @OUT@ -SubPop @WORK@/groups2.txt -MaxDist 100" "$o2.G2.stat.gz" G2.stat.gz

# ---- 9. OutFilterSNP (raw gzip bytes) ----
o=$WORKD/out.filt
"$NEW" -InVCF "$V" -OutStat "$o" -OutFilterSNP -MaxDist 100 >/dev/null 2>&1
emit "filter/vcf" raw "-InVCF @WORK@/sim.vcf.gz -OutStat @OUT@ -OutFilterSNP -MaxDist 100" "$o.vcf.filter.gz" vcf.filter.gz

# ---- 10. InGenotype ----
if [ -f "$GT" ]; then
  o=$WORKD/out.geno
  "$NEW" -InGenotype "$GT" -OutStat "$o" -MaxDist 100 >/dev/null 2>&1
  emit "ingenotype/plain" stat "-InGenotype @WORK@/sim.geno.txt -OutStat @OUT@ -MaxDist 100" "$o.stat.gz" stat.gz
fi

# ---- 11. EHH (raw gzip bytes), on the dense input ----
# sim_ehh.vcf.gz is 10x denser than sim.vcf.gz; -EHH needs more than 8
# non-missing SNPs inside 150 kb of the core site and silently writes an empty
# 20-byte stub otherwise, which would pin nothing at all.
EHH_IN=$WORKD/sim_ehh.vcf.gz
EHH_CORE=$(cat "$WORKD/sim_ehh.core" 2>/dev/null)
if [ -n "$EHH_CORE" ] && [ -f "$EHH_IN" ]; then
  o=$WORKD/out.ehh
  "$NEW" -InVCF "$EHH_IN" -OutStat "$o" -EHH "chr1:$EHH_CORE" >/dev/null 2>&1
  emit "ehh/chr1_core" raw "-InVCF @WORK@/sim_ehh.vcf.gz -OutStat @OUT@ -EHH chr1:$EHH_CORE" "$o.ehh.gz" ehh.gz
  "$NEW" -InVCF "$EHH_IN" -OutStat "$WORKD/out.ehh4" -EHH "chr1:$EHH_CORE" -SubPop "$SUB" >/dev/null 2>&1
  emit "ehh/subpop" raw "-InVCF @WORK@/sim_ehh.vcf.gz -OutStat @OUT@ -EHH chr1:$EHH_CORE -SubPop @WORK@/sub30.txt" "$WORKD/out.ehh4.ehh.gz" ehh.gz
  "$NEW" -InVCF "$EHH_IN" -OutStat "$WORKD/out.ehh5" -EHH "chr1:$EHH_CORE" -OutType 8 >/dev/null 2>&1
  emit "ehh/outtype8" raw "-InVCF @WORK@/sim_ehh.vcf.gz -OutStat @OUT@ -EHH chr1:$EHH_CORE -OutType 8" "$WORKD/out.ehh5.ehh.gz" ehh.gz
fi
# Error paths: an absent contig leaves the same empty stub, and pinning that
# stub is still worth doing -- it is what keeps the failure silent-and-total
# from becoming a crash.
o=$WORKD/out.ehh2
"$NEW" -InVCF "$V" -OutStat "$o" -EHH chrZZ:700000 -MaxDist 100 >/dev/null 2>&1
emit "ehh/absent_chr" raw "-InVCF @WORK@/sim.vcf.gz -OutStat @OUT@ -EHH chrZZ:700000 -MaxDist 100" "$o.ehh.gz" ehh.gz

# ---- 12. phased OutType 0/5 (1,2,3,4,8 already above) ----
# OutTypes 1/6/7/8 all emit the same .stat.gz and differ only in the optional
# pairwise table, so pinning the phased .stat.gz for 6 and 7 as well would add
# rows without adding coverage; their pairwise tables are case 13.
for ot in 0 5; do
  o=$WORKD/out.ph${ot}x
  "$NEW" -InVCF "$VP" -OutStat "$o" -OutType $ot -MaxDist 100 >/dev/null 2>&1
  emit "phased/outtype/$ot" stat "-InVCF @WORK@/sim_phased.vcf.gz -OutStat @OUT@ -OutType $ot -MaxDist 100" "$o.stat.gz" stat.gz
done

# ---- 13. phased pairwise tables ----
for ot in 3 6 7 8; do
  o=$WORKD/out.ph${ot}p
  "$NEW" -InVCF "$VP" -OutStat "$o" -OutType $ot -MaxDist 50 >/dev/null 2>&1
  emit "phased/outtype/$ot/LD" raw "-InVCF @WORK@/sim_phased.vcf.gz -OutStat @OUT@ -OutType $ot -MaxDist 50" "$o.LD.gz" LD.gz
done

# ---- 14. indexed fast path, all OutTypes (-T 4 on bgzip+tabix) ----
for ot in 0 1 2 3 4 5 6 7 8; do
  o=$WORKD/out.fp$ot
  "$NEW" -InVCF "$VG" -OutStat "$o" -OutType $ot -MaxDist 100 -T 4 >/dev/null 2>&1
  emit "fastpath/outtype/$ot" stat "-InVCF @WORK@/sim.bgz.vcf.gz -OutStat @OUT@ -OutType $ot -MaxDist 100 -T 4" "$o.stat.gz" stat.gz
done

# ---- 15. indexed fast path, pairwise tables ----
for ot in 3 6 7 8; do
  o=$WORKD/out.fp${ot}p
  "$NEW" -InVCF "$VG" -OutStat "$o" -OutType $ot -MaxDist 50 -T 4 >/dev/null 2>&1
  emit "fastpath/outtype/$ot/LD" raw "-InVCF @WORK@/sim.bgz.vcf.gz -OutStat @OUT@ -OutType $ot -MaxDist 50 -T 4" "$o.LD.gz" LD.gz
done

# ---- 16. thread extremes on the indexed input ----
for t in 1 16; do
  o=$WORKD/out.t${t}x
  "$NEW" -InVCF "$VG" -OutStat "$o" -T $t -MaxDist 100 >/dev/null 2>&1
  emit "threads/$t" stat "-InVCF @WORK@/sim.bgz.vcf.gz -OutStat @OUT@ -T $t -MaxDist 100" "$o.stat.gz" stat.gz
done

# ---- 17. -L edge cases ----
"$NEW" -InVCF "$VG" -OutStat "$WORKD/out.Lbed" -L "$WORKD/chr1.bed" -T 4 -MaxDist 100 >/dev/null 2>&1
emit "region/bed" stat "-InVCF @WORK@/sim.bgz.vcf.gz -OutStat @OUT@ -L @WORK@/chr1.bed -T 4 -MaxDist 100" "$WORKD/out.Lbed.stat.gz" stat.gz
"$NEW" -InVCF "$VG" -OutStat "$WORKD/out.Ldup" -L chr1,chr1 -T 4 -MaxDist 100 >/dev/null 2>&1
emit "region/dup" stat "-InVCF @WORK@/sim.bgz.vcf.gz -OutStat @OUT@ -L chr1,chr1 -T 4 -MaxDist 100" "$WORKD/out.Ldup.stat.gz" stat.gz
"$NEW" -InVCF "$VG" -OutStat "$WORKD/out.Lnone" -L chrZZ -T 4 -MaxDist 100 >/dev/null 2>&1
emit "region/absent" stat "-InVCF @WORK@/sim.bgz.vcf.gz -OutStat @OUT@ -L chrZZ -T 4 -MaxDist 100" "$WORKD/out.Lnone.stat.gz" stat.gz
if [ -f "$WORKD/chr1only.vcf.gz" ]; then
  "$NEW" -InVCF "$WORKD/chr1only.vcf.gz" -OutStat "$WORKD/out.Lsub" -T 4 -MaxDist 100 >/dev/null 2>&1
  emit "region/chr1_only_subvcf" stat "-InVCF @WORK@/chr1only.vcf.gz -OutStat @OUT@ -T 4 -MaxDist 100" "$WORKD/out.Lsub.stat.gz" stat.gz
fi

# ---- 18. contig without a ##contig line: the fast path must not lose it ----
if [ -f "$WORKD/nocontig.vcf.gz" ]; then
  for t in 1 4; do
    o=$WORKD/out.nc$t
    "$NEW" -InVCF "$WORKD/nocontig.vcf.gz" -OutStat "$o" -T $t -MaxDist 100 >/dev/null 2>&1
    emit "edge/undeclared_contig/T$t" stat "-InVCF @WORK@/nocontig.vcf.gz -OutStat @OUT@ -T $t -MaxDist 100" "$o.stat.gz" stat.gz
  done
fi

# ---- 19. contig name containing ':' -> fast path must decline ----
if [ -f "$WORKD/colon.vcf.gz" ]; then
  for t in 1 4; do
    o=$WORKD/out.cl$t
    "$NEW" -InVCF "$WORKD/colon.vcf.gz" -OutStat "$o" -T $t -MaxDist 100 >/dev/null 2>&1
    emit "edge/colon_contig/T$t" stat "-InVCF @WORK@/colon.vcf.gz -OutStat @OUT@ -T $t -MaxDist 100" "$o.stat.gz" stat.gz
  done
fi

# ---- 20. single-character GT tokens in middle columns ----
if [ -f "$WORKD/shortgt.vcf.gz" ]; then
  "$NEW" -InVCF "$WORKD/shortgt.vcf.gz" -OutStat "$WORKD/out.sg" -MaxDist 100 >/dev/null 2>&1
  emit "edge/short_gt" stat "-InVCF @WORK@/shortgt.vcf.gz -OutStat @OUT@ -MaxDist 100" "$WORKD/out.sg.stat.gz" stat.gz
fi

# ---- 21. SubPop on the indexed fast path, and on phased input ----
"$NEW" -InVCF "$VG" -OutStat "$WORKD/out.subfp" -SubPop "$SUB" -T 4 -MaxDist 100 >/dev/null 2>&1
emit "subpop/indexed" stat "-InVCF @WORK@/sim.bgz.vcf.gz -OutStat @OUT@ -SubPop @WORK@/sub30.txt -T 4 -MaxDist 100" "$WORKD/out.subfp.stat.gz" stat.gz
"$NEW" -InVCF "$VP" -OutStat "$WORKD/out.subph" -SubPop "$SUB" -MaxDist 100 >/dev/null 2>&1
emit "subpop/phased" stat "-InVCF @WORK@/sim_phased.vcf.gz -OutStat @OUT@ -SubPop @WORK@/sub30.txt -MaxDist 100" "$WORKD/out.subph.stat.gz" stat.gz

# ---- 22. repeated -SubPop: two single-column files in one run ----
"$NEW" -InVCF "$V" -OutStat "$WORKD/out.repsub" -SubPop "$WORKD/popA.list" -SubPop "$WORKD/popB.list" -T 1 -MaxDist 100 >/dev/null 2>&1
emit "subpop/repeat/popA" stat "-InVCF @WORK@/sim.vcf.gz -OutStat @OUT@ -SubPop @WORK@/popA.list -SubPop @WORK@/popB.list -T 1 -MaxDist 100" "$WORKD/out.repsub.popA.stat.gz" popA.stat.gz
emit "subpop/repeat/popB" stat "-InVCF @WORK@/sim.vcf.gz -OutStat @OUT@ -SubPop @WORK@/popA.list -SubPop @WORK@/popB.list -T 1 -MaxDist 100" "$WORKD/out.repsub.popB.stat.gz" popB.stat.gz

# ---- 23. InGenotype variants ----
if [ -f "$GT" ]; then
  "$NEW" -InGenotype "$WORKD/sim.geno.plain.gz" -OutStat "$WORKD/out.ggz" -MaxDist 100 >/dev/null 2>&1
  emit "ingenotype/gz" stat "-InGenotype @WORK@/sim.geno.plain.gz -OutStat @OUT@ -MaxDist 100" "$WORKD/out.ggz.stat.gz" stat.gz
  # -SubPop resolves sample names against the header, so this case needs the
  # headered genotype file; against the headerless one it selects nothing and
  # writes an empty .stat.gz.
  "$NEW" -InGenotype "$WORKD/sim.geno.hdr.txt" -OutStat "$WORKD/out.gsub" -SubPop "$WORKD/popA.list" -MaxDist 100 >/dev/null 2>&1
  emit "ingenotype/subpop" stat "-InGenotype @WORK@/sim.geno.hdr.txt -OutStat @OUT@ -SubPop @WORK@/popA.list -MaxDist 100" "$WORKD/out.gsub.stat.gz" stat.gz
  "$NEW" -InGenotype "$GT" -OutStat "$WORKD/out.gmix" -InVCF "$V" -MaxDist 100 >/dev/null 2>&1
  emit "ingenotype/mixed_with_vcf" stat "-InGenotype @WORK@/sim.geno.txt -InVCF @WORK@/sim.vcf.gz -OutStat @OUT@ -MaxDist 100" "$WORKD/out.gmix.stat.gz" stat.gz

  # The PED/MAP route, end to end: make_ci_inputs runs the converter and has
  # already refused to continue if its output differs from the oracle, so these
  # cases freeze what -InGenotype then makes of the converted file. A regression
  # in either the converter or the reader changes the digest.
  if [ -f "$WORKD/sim_numeric.p2g.geno" ]; then
    "$NEW" -InGenotype "$WORKD/sim_numeric.p2g.geno" -OutStat "$WORKD/out.gped" -MaxDist 100 >/dev/null 2>&1
    emit "ingenotype/from_ped_numeric" stat "-InGenotype @WORK@/sim_numeric.p2g.geno -OutStat @OUT@ -MaxDist 100" "$WORKD/out.gped.stat.gz" stat.gz
    "$NEW" -InGenotype "$WORKD/sim_letter.p2g.geno" -OutStat "$WORKD/out.glet" -MaxDist 100 >/dev/null 2>&1
    emit "ingenotype/from_ped_letter" stat "-InGenotype @WORK@/sim_letter.p2g.geno -OutStat @OUT@ -MaxDist 100" "$WORKD/out.glet.stat.gz" stat.gz
    "$NEW" -InGenotype "$WORKD/edge.p2g.geno" -OutStat "$WORKD/out.gedge" -MaxDist 100 >/dev/null 2>&1
    [ -s "$WORKD/out.gedge.stat.gz" ] && \
      emit "ingenotype/from_ped_edge" stat "-InGenotype @WORK@/edge.p2g.geno -OutStat @OUT@ -MaxDist 100" "$WORKD/out.gedge.stat.gz" stat.gz
  fi
  if [ -f "$WORKD/tri.geno" ]; then
    "$NEW" -InGenotype "$WORKD/tri.geno" -OutStat "$WORKD/out.gtri" -SubPop "$WORKD/tri.list" -MaxDist 100 >/dev/null 2>&1
    emit "ingenotype/multiallelic" stat "-InGenotype @WORK@/tri.geno -OutStat @OUT@ -SubPop @WORK@/tri.list -MaxDist 100" "$WORKD/out.gtri.stat.gz" stat.gz
  fi
fi

# ---- 24. -OutStat prefix that already ends in .stat.gz ----
"$NEW" -InVCF "$V" -OutStat "$WORKD/out.prefix.stat.gz" -MaxDist 100 >/dev/null 2>&1
emit "prefix/ends_with_stat_gz" stat "-InVCF @WORK@/sim.vcf.gz -OutStat @OUT@.stat.gz -MaxDist 100" "$WORKD/out.prefix.stat.gz" stat.gz

# ---- 25. contig with a very long name ----------------------------------
# The .LD.gz / .ehh.gz rows start with the chromosome name. The original
# PopLDdecay writes them through a C++ stream, which has no length limit, so
# a long name is exactly the case where this tool has to write the full row:
# the row buffers used to be too small for it, and the length snprintf()
# returns (the length the row WOULD have had) was handed to gzwrite() as the
# write length, so the tail of every row came from beyond the buffer.
if [ -f "$WORKD/longname.vcf" ]; then
  o=$WORKD/out.ln
  "$NEW" -InVCF "$WORKD/longname.vcf" -OutStat "$o" -OutType 1 -MaxDist 100 >/dev/null 2>&1
  emit "edge/long_contig_name/stat" stat "-InVCF @WORK@/longname.vcf -OutStat @OUT@ -OutType 1 -MaxDist 100" "$o.stat.gz" stat.gz
  for ot in 3 8; do
    o=$WORKD/out.ln$ot
    "$NEW" -InVCF "$WORKD/longname.vcf" -OutStat "$o" -OutType $ot -MaxDist 100 >/dev/null 2>&1
    emit "edge/long_contig_name/outtype/$ot" raw "-InVCF @WORK@/longname.vcf -OutStat @OUT@ -OutType $ot -MaxDist 100" "$o.LD.gz" LD.gz
  done
  # the same input on the indexed fast path
  if [ -f "$WORKD/longname.vcf.gz" ]; then
    o=$WORKD/out.lnfp
    "$NEW" -InVCF "$WORKD/longname.vcf.gz" -OutStat "$o" -OutType 8 -MaxDist 100 -T 4 >/dev/null 2>&1
    emit "edge/long_contig_name/fastpath" raw "-InVCF @WORK@/longname.vcf.gz -OutStat @OUT@ -OutType 8 -MaxDist 100 -T 4" "$o.LD.gz" LD.gz
  fi
  # EHH: the table rows start with the chromosome name too
  EHH_CORE=$(cat "$WORKD/longname_ehh.core" 2>/dev/null)
  LONG_NAME=$(cat "$WORKD/longname.name" 2>/dev/null)
  if [ -n "$EHH_CORE" ] && [ -n "$LONG_NAME" ] && [ -f "$WORKD/longname_ehh.vcf.gz" ]; then
    o=$WORKD/out.lnehh
    "$NEW" -InVCF "$WORKD/longname_ehh.vcf.gz" -OutStat "$o" -EHH "$LONG_NAME:$EHH_CORE" >/dev/null 2>&1
    emit "edge/long_contig_name/ehh" raw "-InVCF @WORK@/longname_ehh.vcf.gz -OutStat @OUT@ -EHH $LONG_NAME:$EHH_CORE" "$o.ehh.gz" ehh.gz
  fi
fi

echo
echo "wrote $(grep -vc '^#' "$OUT") fixtures -> $OUT"
