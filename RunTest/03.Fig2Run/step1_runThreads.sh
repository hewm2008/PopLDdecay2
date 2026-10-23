#!/bin/sh
# command-line to run PopLDdecay2 with -T 1/4/8/16/32, five repeats each
#
# The repeats are what the error bars are made of, so do not cut them to one.
# The original PopLDdecay has no -T and is run once per repeat as the baseline.
#
# Each run writes its /bin/time -v record to a .time file; plot_performance.r reads those.
#
# A quick run without downloading chr22 -- uncomment the second line. Absolute
# times are much smaller, the shape of the curves is the same.
IN=../00.GetData/chr22/chr22.vcf.gz
#IN=../00.GetData/base_chr22_50k.vcf.gz

echo Start Time :
date

for r in 1 2 3 4 5; do
for T in 1 4 8 16 32; do
	/bin/time	-v	-p	../../bin/PopLDdecay2	-InVCF	$IN	-OutStat	t$T.r$r	-MaxDist	300	-T	$T	2>	t$T.r$r.time
done
	/bin/time	-v	-p	../../PopLDdecay-3.45/bin/PopLDdecay	-InVCF	$IN	-OutStat	old.r$r	-MaxDist	300	2>	old.r$r.time
	echo "repeat $r done"
done

echo End Time :
date
echo "next: sh step2_makeFigure.sh"
