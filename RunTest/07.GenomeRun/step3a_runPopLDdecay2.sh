#!/bin/sh
# step3: PopLDdecay2 on the whole genome, -T 8
#
# About 33 min and 21 GB peak RSS on an 80-core node. The single-threaded
# baseline takes about 3.7 h and is left commented out below.

G=genome/ALL.autosomes.phase3.vcf.gz

echo Start Time :
date

/bin/time	-v	-p	../../bin/PopLDdecay2	-InVCF	$G	-OutStat	new.t8	-MaxDist	300	-MAF	0.005	-T	8
#/bin/time	-v	-p	../../bin/PopLDdecay2	-InVCF	$G	-OutStat	new.t1	-MaxDist	300	-MAF	0.005	-T	1

echo End Time :
date
