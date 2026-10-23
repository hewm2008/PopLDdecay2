#!/bin/bash
# Run: check PopLDdecay2 against the frozen digests in expected.tsv
#
# 88 synthetic cases. No network, no reference data, no PopLDdecay 3.45 -- this
# is what CI runs. golden_diff.sh runs the same matrix against the real 3.45
# binary instead of frozen digests.
#
# expected.tsv columns: case_id, kind (stat|raw), argv, output_suffix, sha256
#   @WORK@ -> scratch dir, @OUT@ -> output prefix
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)   # this directory
cd "$HERE/../.."                      # repo root

FIX=$HERE/expected.tsv
# The binary under test: bin/ after ./make.sh, build/ when the Makefile
# fallback or a plain cmake build tree is used. NEW= overrides both.
NEW=${NEW:-bin/PopLDdecay2}
[ -x "$NEW" ] || NEW=build/PopLDdecay2
[ -f "$FIX" ] || { echo "missing $FIX"; exit 1; }
[ -x "$NEW" ] || { echo "missing $NEW -- run ./make.sh first"; exit 1; }

WORK=${CI_WORK:-$(mktemp -d)}
mkdir -p "$WORK"
[ -n "${CI_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

# gen_inputs.sh tests [ -f "$BGZIP" ], so a bare command name would silently
# skip the BGZF fixture and fail 26 cases with an empty .stat.gz.
PY=${PYTHON:-python3}
BGZIP=$(command -v "${BGZIP:-bgzip}") || { echo "need bgzip (samtools)"; exit 2; }
TABIX=$(command -v "${TABIX:-tabix}") || { echo "need tabix (samtools)"; exit 2; }
REPO=$(pwd)
export REPO PY BGZIP TABIX
. "$HERE/gen_inputs.sh"
make_ci_inputs "$WORK"

npass=0; nfail=0; nskip=0
while IFS=$'\t' read -r id kind argv suffix want; do
  case "${id:-}" in ''|'#'*) continue ;; esac
  prefix=$WORK/case$(printf '%s' "$id" | tr -c 'A-Za-z0-9' '_')

  # split argv, point -OutStat at our scratch prefix
  set -- ${argv//@WORK@/$WORK}
  out=(); skip=0
  for t in "$@"; do
    [ "$skip" = 1 ] && { out+=("$prefix"); skip=0; continue; }
    [ "$t" = "-OutStat" ] && { out+=("-OutStat"); skip=1; continue; }
    out+=("${t//@OUT@/$prefix}")
  done
  "$NEW" "${out[@]}" >/dev/null 2>&1

  produced=$prefix.$suffix
  [ -f "$produced" ] || {
    printf 'SKIP  %-34s (%s not produced)\n' "$id" "$suffix"
    nskip=$((nskip+1)); continue
  }
  if [ "$kind" = stat ]; then
    got=$(zcat "$produced" 2>/dev/null | sha256sum | cut -d' ' -f1)
  else
    got=$(sha256sum "$produced" 2>/dev/null | cut -d' ' -f1)
  fi

  if [ "$got" = "$want" ]; then
    printf 'PASS  %-34s %s\n' "$id" "${want:0:16}"; npass=$((npass+1))
  else
    printf 'FAIL  %-34s\n        want %s\n        got  %s\n' "$id" "$want" "$got"
    nfail=$((nfail+1))
  fi
done < "$FIX"

# ---------------------------------------------------------------------------
# Behaviour cases: exit status and messages, not digests.
#
# The digest table above pins what a successful run writes. These pin what a
# run that cannot succeed does: an invalid value must be reported instead of
# being silently reinterpreted (atof("abc") is 0, which switched the MAF
# filter off; atoi("abc") is 0, which turned -OutType into a header-only
# table), and a failed write must not leave a truncated file behind.
# ---------------------------------------------------------------------------
ncase=0   # behaviour cases run below; the doc check counts them

# expect_fail ID EXPECTED_STDERR_SUBSTRING -- ARGS...
expect_fail() {
  local id=$1 want_msg=$2; shift 2
  local err="$WORK/err.$(printf '%s' "$id" | tr -c 'A-Za-z0-9' '_')"
  ncase=$((ncase+1))
  "$NEW" "$@" >/dev/null 2>"$err"
  local rc=$?
  if [ "$rc" -eq 0 ]; then
    printf 'FAIL  %-34s (exited 0, expected a failure)\n' "$id"; nfail=$((nfail+1)); return
  fi
  if [ -n "$want_msg" ] && ! grep -qF -- "$want_msg" "$err"; then
    printf 'FAIL  %-34s (message %s not reported)\n' "$id" "\"$want_msg\""; nfail=$((nfail+1)); return
  fi
  printf 'PASS  %-34s (exit %d)\n' "$id" "$rc"; npass=$((npass+1))
}

# expect_ok ID -- ARGS...
expect_ok() {
  local id=$1; shift
  local err="$WORK/err.$(printf '%s' "$id" | tr -c 'A-Za-z0-9' '_')"
  ncase=$((ncase+1))
  "$NEW" "$@" >/dev/null 2>"$err"
  local rc=$?
  if [ "$rc" -ne 0 ]; then
    printf 'FAIL  %-34s (exit %d, expected success)\n' "$id" "$rc"; nfail=$((nfail+1)); return
  fi
  printf 'PASS  %-34s (exit 0)\n' "$id"; npass=$((npass+1))
}

SIM="$WORK/sim.vcf.gz"
OUTP="$WORK/behave"

# Invalid numbers must be refused, not reinterpreted.
expect_fail cli/maf_not_a_number   "-MAF should be a number"   -InVCF "$SIM" -OutStat "$OUTP" -MAF abc
expect_fail cli/maf_nan             "-MAF should be a number"   -InVCF "$SIM" -OutStat "$OUTP" -MAF nan
expect_fail cli/maf_out_of_range    "-MAF should be a number"   -InVCF "$SIM" -OutStat "$OUTP" -MAF 2
expect_fail cli/maf_negative        "-MAF should be a number"   -InVCF "$SIM" -OutStat "$OUTP" -MAF -0.5
expect_fail cli/het_out_of_range    "-Het should be a number"   -InVCF "$SIM" -OutStat "$OUTP" -Het 1.5
expect_fail cli/miss_not_a_number   "-Miss should be a number"  -InVCF "$SIM" -OutStat "$OUTP" -Miss xyz
expect_fail cli/maxdist_not_int     "-MaxDist should be an integer" -InVCF "$SIM" -OutStat "$OUTP" -MaxDist 1e9
expect_fail cli/maxdist_trailing    "-MaxDist should be an integer" -InVCF "$SIM" -OutStat "$OUTP" -MaxDist 100kb
expect_fail cli/outtype_not_int     "-OutType should be an integer" -InVCF "$SIM" -OutStat "$OUTP" -OutType abc
expect_fail cli/outtype_out_of_range "-OutType should be an integer" -InVCF "$SIM" -OutStat "$OUTP" -OutType 9
expect_fail cli/threads_not_int     "-T should be an integer"   -InVCF "$SIM" -OutStat "$OUTP" -T 8x
expect_fail cli/method_not_int      "-Method should be an integer" -InVCF "$SIM" -OutStat "$OUTP" -Method two
expect_fail cli/missing_value       "Lack Argument"             -InVCF "$SIM" -OutStat -MaxDist 100
expect_fail cli/ehh_not_a_number    "the site must be a positive whole number" -InVCF "$SIM" -OutStat "$OUTP" -EHH chr1:abc

# The boundary values that the author's own scripts use must keep working:
# -MAF 0 / -Miss 1 / -Het 1 open every filter up (rebuttal/03.correctness).
expect_ok cli/boundary_filters -InVCF "$SIM" -OutStat "$OUTP.b" -MAF 0 -Miss 1 -Het 1 -MaxDist 100

# An output that cannot be written must fail, and must not leave a file that
# looks like a result.
ncase=$((ncase+1))
mkdir -p "$WORK/nowrite"
chmod a-w "$WORK/nowrite" 2>/dev/null
if [ -w "$WORK/nowrite" ]; then
  printf 'SKIP  %-34s (the scratch directory is writable anyway)\n' "io/unwritable_output"
  nskip=$((nskip+1))
else
  if "$NEW" -InVCF "$SIM" -OutStat "$WORK/nowrite/x" -MaxDist 10 >/dev/null 2>&1; then
    printf 'FAIL  %-34s (exited 0)\n' "io/unwritable_output"; nfail=$((nfail+1))
  else
    printf 'PASS  %-34s (exit non-zero)\n' "io/unwritable_output"; npass=$((npass+1))
  fi
  chmod u+w "$WORK/nowrite" 2>/dev/null
fi

# The EHH path keeps the original's exit status: a run that cannot find the
# chromosome still exits 0 (the original does the same), so batch jobs that
# compare exit codes against 3.45 keep working. Pinned here so the
# compatibility constraint cannot be broken by accident.
ncase=$((ncase+1))
"$NEW" -InVCF "$SIM" -OutStat "$OUTP.ehh" -EHH chrZZ:700000 >/dev/null 2>&1
if [ $? -eq 0 ]; then
  printf 'PASS  %-34s (exit 0, as in 3.45)\n' "ehh/absent_chrom_exit_0"; npass=$((npass+1))
else
  printf 'FAIL  %-34s (exit status changed; 3.45 exits 0 here)\n' "ehh/absent_chrom_exit_0"; nfail=$((nfail+1))
fi

# ---------------------------------------------------------------------------
# Documentation cases: the two READMEs must not drift apart.
#
# The English and the Chinese README are written by hand and separately, so
# they drift: an option was added to the English table and not to the Chinese
# one (-InPlink was missing for a whole release). These checks catch that.
# ---------------------------------------------------------------------------

# The option names listed in a README's parameter block (the tool's own -help
# output), one per line.
doc_options() {   # readme
  sed -n '/Usage: PopLDdecay2/,/Show more help/p' "$1" \
    | grep -oE '^[[:space:]]+-[A-Za-z]+' | tr -d ' \t' | sort -u
}

ncase=$((ncase+1))
en_opts=$(doc_options README.md)
cn_opts=$(doc_options README.cn.md)
if [ "$en_opts" = "$cn_opts" ]; then
  printf 'PASS  %-34s (%s options)\n' "docs/option_tables_match" "$(echo "$en_opts" | wc -l)"
  npass=$((npass+1))
else
  printf 'FAIL  %-34s\n' "docs/option_tables_match"
  diff <(echo "$en_opts") <(echo "$cn_opts") | head -10
  nfail=$((nfail+1))
fi

# Every relative link in either README must resolve.
ncase=$((ncase+1))
doc_links_ok=1
for readme in README.md README.cn.md; do
  for link in $(grep -oE '\]\([^)]+\)' "$readme" | sed 's/^](//; s/)$//' | sort -u); do
    case "$link" in
      http*|mailto:*|''|'#'*) continue ;;   # '#...' is an in-document anchor
    esac
    target="${link%%#*}"
    [ -e "$target" ] || { printf '        %s: broken link -> %s\n' "$readme" "$link"; doc_links_ok=0; }
  done
done
if [ "$doc_links_ok" = 1 ]; then
  printf 'PASS  %-34s (all relative links resolve)\n' "docs/links_resolve"
  npass=$((npass+1))
else
  printf 'FAIL  %-34s\n' "docs/links_resolve"; nfail=$((nfail+1))
fi

# The number of cases the READMEs advertise must be the number this script
# actually runs: the digest cases in expected.tsv plus the behaviour cases
# above.
ncase=$((ncase+1))   # this case counts itself
n_digest=$(grep -vc '^#' "$FIX")
n_total=$((n_digest + ncase))
docs_ok=1
for readme in README.md README.cn.md; do
  # A digit boundary, so "105" does not match inside "1050".
  grep -qE "(^|[^0-9])$n_total([^0-9]|$)" "$readme" \
    || { printf '        %s: does not mention %s cases\n' "$readme" "$n_total"; docs_ok=0; }
done
if [ "$docs_ok" = 1 ]; then
  printf 'PASS  %-34s (both READMEs say %s cases)\n' "docs/case_count_matches" "$n_total"
  npass=$((npass+1))
else
  printf 'FAIL  %-34s (expected %s cases)\n' "docs/case_count_matches" "$n_total"; nfail=$((nfail+1))
fi

echo
echo "passed: $npass   failed: $nfail   skipped: $nskip"
[ "$nfail" -eq 0 ] || exit 1
