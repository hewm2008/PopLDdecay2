#!/bin/sh
# step3b: the pairwise table of the three programs, side by side
#
# The three programs write the same pairs in the same order, so the first rows
# of the three tables line up and can be read against each other:
#
#   PLINK 2.0   output/plink/plink2_<TAG>.vcor   -> chrom, pos1, pos2, PHASED_R2
#   Haploview   output/haploview/<TAG>.LD        -> L1, L2, r^2
#   PopLDdecay2 output/poplddecay2/<TAG>.LD.gz   -> chr, Site1, Site2, r^2
#
# This is a look, not a comparison: the numbers are not compared against each
# other here, and not against the ground truth. step3_compare.sh does that, and
# writes output/compare/<TAG>/report.txt.
#
#   sh step3b_showPairs.sh

#echo Start Time :
#date

TAG=${TAG:-sub100_maf0.1}

echo Plink2
head  output/plink/plink2_$TAG.vcor  | cut -f 1,2,4,5
echo Haploview4.2
head  output/haploview/$TAG.LD       | cut -f 1,2,5
echo PopLDdecay2
zcat  output/poplddecay2/$TAG.LD.gz  | head | cut -f 1-3,6

#echo End Time :
#date
