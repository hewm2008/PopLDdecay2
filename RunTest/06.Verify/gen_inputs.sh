#!/bin/bash
# Build the synthetic inputs used by step1_checkOutput.sh and make_fixtures.sh.
#
#   . gen_inputs.sh ; make_ci_inputs "$WORK"
#
# Both scripts source this file, so the inputs a digest was frozen from are
# byte-for-byte the inputs the check rebuilds. Everything is deterministic:
# gen_sim_vcf.py is seeded and every derived file follows from the ones before.
#
# Requires BGZIP, TABIX, PY and the gen_sim_*.py next to this file.
make_ci_inputs() {   # workdir
  local W=$1
  mkdir -p "$W"
  local PY=${PY:-python3}
  local GEN="${GEN:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
  # A 200-character chromosome name: long enough to overflow every row
  # buffer in the .LD.gz / .ehh.gz writers (the shortest is 128 bytes, and
  # the numeric part of a row is ~35 bytes).
  local LONG_CHROM_NAME="chr$(printf 'L%.0s' {1..197})"
  # Written out so the other scripts can use the same name in a command line.
  [ -f "$W/longname.name" ] || printf '%s\n' "$LONG_CHROM_NAME" > "$W/longname.name"

  # ---- primary inputs -------------------------------------------------
  [ -f "$W/sim.vcf" ] || "$PY" "$GEN/gen_sim_vcf.py" -o "$W/sim.vcf" \
                             --nsamples 60 --sites 2000 --chroms 4 >/dev/null 2>&1
  [ -f "$W/sim.vcf.gz" ] || gzip -c "$W/sim.vcf" > "$W/sim.vcf.gz"
  [ -f "$W/sim_phased.vcf.gz" ] || "$PY" "$GEN/gen_sim_vcf.py" -o "$W/sim_phased.vcf.gz" \
                             --gzip --phased --nsamples 60 --sites 2000 --chroms 4 >/dev/null 2>&1
  if [ -f "$BGZIP" ] && [ -f "$W/sim.vcf" ]; then
    [ -f "$W/sim.bgz.vcf.gz" ] || "$BGZIP" -f -c "$W/sim.vcf" > "$W/sim.bgz.vcf.gz"
    [ -f "$W/sim.bgz.vcf.gz.tbi" ] && [ -f "$TABIX" ] && \
      "$TABIX" -f -p vcf "$W/sim.bgz.vcf.gz" 2>/dev/null
  fi
  [ -f "$W/sim.geno.txt" ] || "$PY" "$GEN/gen_sim_genotype.py" -i "$W/sim.vcf" \
                             -o "$W/sim.geno.txt" --nsamples 60 >/dev/null 2>&1
  [ -f "$W/sim.geno.plain.gz" ] || gzip -c "$W/sim.geno.txt" > "$W/sim.geno.plain.gz"
  # -SubPop against -InGenotype needs the sample-name header, which the plain
  # generator does not emit.
  [ -f "$W/sim.geno.hdr.txt" ] || "$PY" "$GEN/gen_sim_genotype.py" -i "$W/sim.vcf" \
                             -o "$W/sim.geno.hdr.txt" --nsamples 60 --with-header >/dev/null 2>&1

  # ---- PED/MAP for bin/mis/plink2genotype.pl ---------------------------
  # The converter turns PED/MAP into the -InGenotype format, and until now it
  # had no test anywhere: the golden suite and the CI gate both exercised
  # -InGenotype only through gen_sim_genotype.py, which bypasses it entirely.
  # gen_sim_ped.py writes the PED/MAP and the .geno they must convert to, so the
  # converter can be checked against an oracle rather than against itself.
  #
  # Running it here, during input generation, means a broken converter fails
  # both the fixture maker and CI, loudly, instead of silently producing a
  # wrong input that the -InGenotype cases below would then bless.
  for enc in numeric letter; do
    [ -f "$W/sim_$enc.ped" ] || "$PY" "$GEN/gen_sim_ped.py" -i "$W/sim.vcf" \
        -o "$W/sim_$enc" --encoding "$enc" >/dev/null 2>&1
    [ -f "$W/sim_$enc.p2g.geno" ] || {
      perl "$REPO/bin/mis/plink2genotype.pl" \
           -inPED "$W/sim_$enc.ped" -inMAP "$W/sim_$enc.map" \
           -outGenotype "$W/sim_$enc.p2g.geno" >/dev/null 2>&1
      if ! cmp -s "$W/sim_$enc.p2g.geno" "$W/sim_$enc.geno"; then
        echo "FATAL: bin/mis/plink2genotype.pl ($enc PED) does not reproduce" >&2
        echo "       the .geno written by gen_sim_ped.py" >&2
        echo "       diff: diff $W/sim_$enc.geno $W/sim_$enc.p2g.geno" >&2
        return 1
      fi
    }
  done
  [ -f "$W/edge.ped" ] || "$PY" "$GEN/gen_sim_ped.py" --edge -o "$W/edge" \
      --encoding letter >/dev/null 2>&1
  [ -f "$W/edge.p2g.geno" ] || {
    perl "$REPO/bin/mis/plink2genotype.pl" \
         -inPED "$W/edge.ped" -inMAP "$W/edge.map" \
         -outGenotype "$W/edge.p2g.geno" >/dev/null 2>&1
    if ! cmp -s "$W/edge.p2g.geno" "$W/edge.geno"; then
      echo "FATAL: bin/mis/plink2genotype.pl (edge fixture) does not reproduce" >&2
      echo "       gen_sim_ped.py --edge output" >&2
      return 1
    fi
  }

  # -EHH needs more than 8 non-missing SNPs within 150 kb of the core site.
  # sim.vcf (2,000 sites over 4 chromosomes) is far too sparse for that, so EHH
  # gets its own dense input and a core site picked out of it.
  if [ ! -f "$W/sim_ehh.vcf.gz" ]; then
    "$PY" "$GEN/gen_sim_vcf.py" -o "$W/sim_ehh.raw.vcf.gz" --gzip --phased \
          --nsamples 60 --sites 20000 --chroms 2 >/dev/null 2>&1
    # The generator deliberately emits missing calls (".", "0/."), and -EHH
    # silently writes a 20-byte stub when fewer than 8 SNPs in the window are
    # fully called. Complete them here -- GT is the only FORMAT field, so
    # turning each missing allele into "0" is a well-defined rewrite.
    { zcat "$W/sim_ehh.raw.vcf.gz" | awk 'BEGIN{FS=OFS="\t"}
        /^#/ {print; next}
        { for (i=10; i<=NF; i++) gsub(/\./, "0", $i); print }'; } \
      | "$BGZIP" -c > "$W/sim_ehh.vcf.gz"
    [ -f "$TABIX" ] && "$TABIX" -f -p vcf "$W/sim_ehh.vcf.gz" 2>/dev/null
  fi
  [ -f "$W/sim_ehh.core" ] || zcat "$W/sim_ehh.vcf.gz" 2>/dev/null \
      | awk -F'\t' '$1=="chr1" && !/^#/{n++; if(n==200) print $2}' > "$W/sim_ehh.core"

  # The generator keeps every sample on the single #CHROM line.
  samples() { zcat "$W/sim.vcf.gz" 2>/dev/null | grep -m1 '^#CHROM' | cut -f10- | tr '\t' '\n'; }

  # ---- population lists ----------------------------------------------
  [ -f "$W/sub30.txt" ]     || samples | head -30 > "$W/sub30.txt"
  [ -f "$W/groups2.txt" ]   || samples | head -40 | awk '{print $1"\t"((NR<=20)?"G1":"G2")}' > "$W/groups2.txt"
  [ -f "$W/popA.list" ]     || samples | sed -n '1,15p'  > "$W/popA.list"
  [ -f "$W/popB.list" ]     || samples | sed -n '16,30p' > "$W/popB.list"

  # ---- -L BED (0-based, half-open) -----------------------------------
  [ -f "$W/chr1.bed" ] || printf 'chr1\t0\t400000\n' > "$W/chr1.bed"

  # ---- chr1-only sub-VCF: -L chr1 on the full file must equal this ----
  if [ ! -f "$W/chr1only.vcf.gz" ] && [ -f "$BGZIP" ]; then
    { zcat "$W/sim.bgz.vcf.gz" | awk 'BEGIN{FS=OFS="\t"} /^#/{print;next} $1=="chr1"{print}' \
      | "$BGZIP" -c; } > "$W/chr1only.vcf.gz"
    [ -f "$TABIX" ] && "$TABIX" -f -p vcf "$W/chr1only.vcf.gz" 2>/dev/null
  fi

  # ---- contig with NO ##contig line (fastpath must not drop its data) --
  if [ ! -f "$W/nocontig.vcf.gz" ] && [ -f "$BGZIP" ]; then
    { zcat "$W/sim.bgz.vcf.gz" | awk 'BEGIN{FS=OFS="\t"}
        /^##contig/{next} {print}' | "$BGZIP" -c; } > "$W/nocontig.vcf.gz"
    [ -f "$TABIX" ] && "$TABIX" -f -p vcf "$W/nocontig.vcf.gz" 2>/dev/null
  fi

  # ---- contig name containing ':' (tabix query is unsafe -> decline) ---
  if [ ! -f "$W/colon.vcf.gz" ] && [ -f "$BGZIP" ]; then
    { zcat "$W/sim.bgz.vcf.gz" | awk 'BEGIN{FS=OFS="\t"}
        /^#CHROM/{sub(/\tchr1\t/,"\tchr:1\t"); print; next}
        /^#/{next}
        $1=="chr1"{$1="chr:1"; print}' | "$BGZIP" -c; } > "$W/colon.vcf.gz"
    [ -f "$TABIX" ] && "$TABIX" -f -p vcf "$W/colon.vcf.gz" 2>/dev/null
  fi

  # ---- single-character GT tokens in the MIDDLE sample columns --------
  # A haploid "0" is one character wide, so a naive third-character read picks
  # up the next field. Columns 14-18 become one- or two-character tokens while
  # the last column is left untouched, so the case does not depend on where
  # the line happens to end.
  [ -f "$W/shortgt.vcf.gz" ] || {
    zcat "$W/sim.vcf.gz" 2>/dev/null | awk 'BEGIN{FS=OFS="\t"}
      /^#/ {print; next}
      { for (i=14; i<=18; i++) $i = ($i == "0/0" || $i == "0|0") ? "0" : "1|0"
        print }' | gzip -c > "$W/shortgt.vcf.gz"
  }

  # ---- contig with a very long name ------------------------------------
  # The .LD.gz / .ehh.gz row buffers are sized for a short chromosome name.
  # A long one used to overflow them: snprintf() returns the length the row
  # WOULD have had, and that length was handed to gzwrite() as the write
  # length, so the tail of every row was read from beyond the buffer. The
  # original PopLDdecay wrote these rows through a C++ stream, which has no
  # such limit, so the long name is the case where the two must agree.
  [ -f "$W/longname.vcf" ] || "$PY" "$GEN/gen_sim_vcf.py" -o "$W/longname.vcf" \
      --chrom-name "$LONG_CHROM_NAME" --nsamples 60 --sites 400 --chroms 1 >/dev/null 2>&1
  if [ ! -f "$W/longname.vcf.gz" ] && [ -f "$BGZIP" ]; then
    "$BGZIP" -f -c "$W/longname.vcf" > "$W/longname.vcf.gz"
    [ -f "$TABIX" ] && "$TABIX" -f -p vcf "$W/longname.vcf.gz" 2>/dev/null
  fi

  # ---- the same dense EHH input, under the long chromosome name ---------
  # The EHH table starts every row with the chromosome name, so it needs the
  # long-name case as much as the .LD.gz table does.
  if [ ! -f "$W/longname_ehh.vcf.gz" ] && [ -f "$BGZIP" ] && [ -f "$W/sim_ehh.vcf.gz" ]; then
    { zcat "$W/sim_ehh.vcf.gz" | awk -v nm="$LONG_CHROM_NAME" 'BEGIN{FS=OFS="\t"}
        /^##contig=<ID=chr1,/{sub(/ID=chr1,/, "ID=" nm ","); print; next}
        /^#/{print; next}
        $1=="chr1"{$1=nm; print}
        {print}' | "$BGZIP" -c; } > "$W/longname_ehh.vcf.gz"
    [ -f "$TABIX" ] && "$TABIX" -f -p vcf "$W/longname_ehh.vcf.gz" 2>/dev/null
  fi
  [ -f "$W/longname_ehh.core" ] || zcat "$W/longname_ehh.vcf.gz" 2>/dev/null \
      | awk -F'\t' -v nm="$LONG_CHROM_NAME" \
        'BEGIN{k=0} $1==nm && !/^#/{k++; if(k==200) print $2}' > "$W/longname_ehh.core"

  # ---- multi-allelic + singleton genotype input ------------------------
  [ -f "$W/tri.geno" ] || printf '#CHROM\tPOS\tS1 S2 S3 S4\nchr1\t100\tA A A A\nchr1\t200\tA C G A\nchr1\t300\tA C A C\n' > "$W/tri.geno"
  [ -f "$W/tri.list" ] || printf 'S1 S2 S3 S4\n' > "$W/tri.list"
}