#!/bin/sh
# command-line to run PopLDdecay2, one run per population
#
# The baseline of step1a: the same tool, the same five populations, but five
# separate runs with a single-column sample list instead of one pass with the
# two-column sample.group. The outputs must be byte-identical to the ones of
# step1a_runPopLDdecay2Multi.sh (same sites, same pairs) -- that is what
# step3_checkOutput.sh checks -- and the wall times say what the one-pass
# design saves.
#
# -T 8 is the same as in step1a, so the two are comparable.

IN=../00.GetData/chr22/chr22.vcf.gz
#IN=../00.GetData/base_chr22_50k.vcf.gz

for p in AFR AMR EAS EUR SAS; do
	awk -v p=$p '$2==p {print $1}' sample.group > $p.list
done
wc -l *.list

echo Start Time :
date

/bin/time	-v	-p	../../bin/PopLDdecay2	-InVCF	$IN	-OutStat	AFR.NewS.stat	-SubPop	AFR.list	-MaxDist	300	-T	8	2>	AFR.NewS.time
/bin/time	-v	-p	../../bin/PopLDdecay2	-InVCF	$IN	-OutStat	AMR.NewS.stat	-SubPop	AMR.list	-MaxDist	300	-T	8	2>	AMR.NewS.time
/bin/time	-v	-p	../../bin/PopLDdecay2	-InVCF	$IN	-OutStat	EAS.NewS.stat	-SubPop	EAS.list	-MaxDist	300	-T	8	2>	EAS.NewS.time
/bin/time	-v	-p	../../bin/PopLDdecay2	-InVCF	$IN	-OutStat	EUR.NewS.stat	-SubPop	EUR.list	-MaxDist	300	-T	8	2>	EUR.NewS.time
/bin/time	-v	-p	../../bin/PopLDdecay2	-InVCF	$IN	-OutStat	SAS.NewS.stat	-SubPop	SAS.list	-MaxDist	300	-T	8	2>	SAS.NewS.time

echo End Time :
date
echo "next: sh step3_checkOutput.sh"
