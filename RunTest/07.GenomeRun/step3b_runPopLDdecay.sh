#!/bin/sh
# step4: the original PopLDdecay on the whole genome
#
# A full run takes about 61 h with a 75 GB peak RSS on an 80-core node, so this
# is the slow step by far. The 60-second smoke test below is commented out: it
# is enough to confirm that the process starts, reads, and keeps working,
# without occupying the machine for days. Exit code 124 from it means the
# timeout fired, which is the expected result, not a failure.

G=genome/ALL.autosomes.phase3.vcf.gz

echo Start Time :
date

# smoke test (60 s):
#/bin/time	-v	-p	timeout 60 ../../PopLDdecay-3.45/bin/PopLDdecay	-InVCF	$G	-OutStat	old	-MaxDist	300	-MAF	0.005

# full run (about 61 h):
/bin/time	-v	-p	../../PopLDdecay-3.45/bin/PopLDdecay	-InVCF	$G	-OutStat	old	-MaxDist	300	-MAF	0.005

echo End Time :
date
