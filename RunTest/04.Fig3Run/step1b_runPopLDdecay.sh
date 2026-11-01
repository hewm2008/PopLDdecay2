#!/bin/sh
# command-line to run the original PopLDdecay, one run per population
#
# It accepts only a single-column sample list, so sample.group is split into five
# lists and the tool runs five times, sequentially. step1a_runPopLDdecay2Multi.sh
# does the same work in ONE pass; step3_checkOutput.sh compares the two output
# sets and the two wall times.
#
# PopLDdecay is not redistributable, so it is not in this package. Expected next
# in the package root: PopLDdecay-3.45/bin/PopLDdecay

IN=../00.GetData/chr22/chr22.vcf.gz
#IN=../00.GetData/base_chr22_50k.vcf.gz

for p in AFR AMR EAS EUR SAS; do
	awk -v p=$p '$2==p {print $1}' sample.group > $p.list
done
wc -l *.list

echo Start Time :
date

/bin/time	-v	-p	../../PopLDdecay-3.45/bin/PopLDdecay	-InVCF	$IN	-OutStat	AFR.old.stat	-SubPop	AFR.list	-MaxDist	300	2>	AFR.old.time
/bin/time	-v	-p	../../PopLDdecay-3.45/bin/PopLDdecay	-InVCF	$IN	-OutStat	AMR.old.stat	-SubPop	AMR.list	-MaxDist	300	2>	AMR.old.time
/bin/time	-v	-p	../../PopLDdecay-3.45/bin/PopLDdecay	-InVCF	$IN	-OutStat	EAS.old.stat	-SubPop	EAS.list	-MaxDist	300	2>	EAS.old.time
/bin/time	-v	-p	../../PopLDdecay-3.45/bin/PopLDdecay	-InVCF	$IN	-OutStat	EUR.old.stat	-SubPop	EUR.list	-MaxDist	300	2>	EUR.old.time
/bin/time	-v	-p	../../PopLDdecay-3.45/bin/PopLDdecay	-InVCF	$IN	-OutStat	SAS.old.stat	-SubPop	SAS.list	-MaxDist	300	2>	SAS.old.time

echo End Time :
date
echo "next: sh step3_checkOutput.sh"
