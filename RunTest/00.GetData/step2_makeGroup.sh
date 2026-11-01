#!/bin/sh
# step2: make sample.group (sample <TAB> super-population) from the 1000G panel
#
# sample.group is the input of -SubPop: with this two-column file one PopLDdecay2
# run outputs one .stat.gz per population. The original PopLDdecay takes only a
# single-column list, so it needs one run per population.
#
# The result is already shipped as 04.Fig3Run/sample.group; run this only if you
# want to rebuild it from your own panel.

PANEL=${1:-chr22/panel.txt}

awk 'BEGIN{OFS="\t"} !/^#/ && NF>=3 {print $1, $3}' $PANEL > sample.group

echo "wrote sample.group : $(wc -l < sample.group) samples"
cut -f2 sample.group | sort | uniq -c
