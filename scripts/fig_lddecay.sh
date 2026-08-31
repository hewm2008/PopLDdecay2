#!/bin/bash
# Generate all paper figures into paper/Fig (PDF+PNG, vector) from benchmark data.
#   Fig1_workflow       : graphviz dot (work flow)
#   Fig2_chr22_lddecay  : whole-panel chr22 LD decay (Plot_OnePop.pl)
#   Fig3_AFR..SAS       : per super-population LD decay (Plot_OnePop.pl)
#   Fig3_subpop         : all five super-populations overlaid (Plot_MultiPop.pl)
#   Fig4_{scaling,speedup,mem} : performance figures (scripts/fig_perf.R)
# Usage: bash scripts/fig_lddecay.sh [workdir] [outdir]
#   workdir defaults to Test//paper_chr22 (whole panel) or paper_subpop (groups).
#   Whole-panel: uses old.345/new.t1 .stat.gz; SubPop: new.<GROUP>.stat.gz.
set -uo pipefail
cd "$(dirname "$0")/.."
PLOT=PopLDdecay-3.45/bin/Plot_OnePop.pl
MPLOT=PopLDdecay-3.45/bin/Plot_MultiPop.pl
OUT=${2:-paper/Fig}
DOTSRC=${3:-paper/Fig/Fig1_workflow.dot}

mkdir -p "$OUT"

# ---- Figure 1: workflow (graphviz) ----
if [ -f "$DOTSRC" ]; then
  dot -Tpdf "$DOTSRC" -o "$OUT/Fig1_workflow.pdf" && \
  dot -Tpng -Gdpi=300 "$DOTSRC" -o "$OUT/Fig1_workflow.png" && \
    echo "fig1: $OUT/Fig1_workflow.{pdf,png}" || echo "fig1 FAILED"
fi

# ---- Figure 2: whole chr22 LD decay (3.45 vs PopLDdecay2 overlay is via MutiPop) ----
# Whole panel single curve from new.t1 (byte-identical to old.345 anyway).
W=Test//paper_chr22
if [ -f "$W/new.t1.stat.gz" ]; then
  perl "$PLOT" -inFile "$W/new.t1.stat.gz" -output "$OUT/Fig2_chr22_lddecay" -maxX 300 2>/dev/null \
    && echo "fig2: $OUT/Fig2_chr22_lddecay.{pdf,png}" || echo "fig2 FAILED"
fi

# ---- Figure 3: multi-population LD decay (AFR/AMR/EAS/EUR/SAS) ----
S=Test//paper_subpop
LIST=Test//fig_subpop.list
: > "$LIST"
if [ -d "$S" ]; then
  for g in AFR AMR EAS EUR SAS; do
    if [ -f "$S/new.$g.stat.gz" ]; then
      perl "$PLOT" -inFile "$S/new.$g.stat.gz" -output "$OUT/Fig3_$g" -maxX 300 2>/dev/null \
        && echo "fig3/$g done" || echo "fig3/$g FAILED"
      echo "$S/new.$g.stat.gz $g" >> "$LIST"
    fi
  done
  if [ -s "$LIST" ]; then
    perl "$MPLOT" -inList "$LIST" -output "$OUT/Fig3_subpop" 2>/dev/null \
      && echo "fig3/subpop: $OUT/Fig3_subpop.{pdf,png}" || echo "fig3/subpop FAILED"
  fi
fi

# ---- Figure 4: performance (thread scaling / speedup / memory) ----
if command -v Rscript >/dev/null 2>&1; then
  Rscript scripts/fig_perf.R Test//paper_chr22 "$OUT/Fig4" >/dev/null 2>&1 \
    && echo "fig4: $OUT/Fig4_{scaling,speedup,mem}.{pdf,png}" || echo "fig4 FAILED"
else
  echo "fig4 SKIPPED (Rscript not available)"
fi

echo "done -> $OUT"
