#!/bin/sh
# command-line to run PopLDdecay2 once for five super-populations
#
# sample.group is a two-column file (sample <TAB> super-population), so ONE
# indexed pass writes one .stat.gz per population. The original PopLDdecay
# accepts only a single-column list and needs five runs -- see step1b_runPopLDdecay.sh.
#
# A quick run without downloading chr22 -- uncomment the second line.
IN=../00.GetData/chr22/chr22.vcf.gz
#IN=../00.GetData/base_chr22_50k.vcf.gz

echo Start Time :
date

/bin/time	-v	-p	../../bin/PopLDdecay2	-InVCF	$IN	-OutStat	mp	-SubPop	sample.group	-MaxDist	300	-T	8	2>	mp.time

ls -1 mp.*.stat.gz

echo End Time :
date
echo "next: step1b_runPopLDdecay.sh is the same work with the old tool; step2_plotCurves.sh draws Figure 3"
