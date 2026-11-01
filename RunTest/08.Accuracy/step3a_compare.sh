#!/bin/sh
# step3a: compare the pairwise R^2 (and D', LOD, CI) of the three programs
#
#   tools/ld_truth.py     the ground truth: count the phased haplotypes of the
#                         VCF directly and apply the textbook formulas
#   tools/compare_ld.py   join the three tables on the site pair and report
#
# Everything is written to output/compare/<TAG>/:
#
#   report.txt        the summary that is worth reading
#   r2_pairs.tsv      THE TABLE: one row per site pair, the R^2 of every program
#   r2_scatter.png    R^2 of the programs against each other and against the truth
#   pairwise.tsv      the same table with D', LOD and CI added
#   curve.tsv         Mean_r^2 per 1 bp distance bin
#   curve_10kb.tsv    Mean_r^2 per 10 kb distance bin
#   truth_ld.tsv      the ground truth table
#
#   sh step3a_compare.sh

echo Start Time :
date

TAG=${TAG:-sub100_maf0.1}

IN=input/$TAG.vcf.gz
DIR=output/compare/$TAG

mkdir -p $DIR

python3 tools/ld_truth.py $IN $DIR/truth_ld.tsv 300000
python3 tools/compare_ld.py $TAG

echo End Time :
date
echo "next: sh step4a_runUnphased.sh"
