#!/bin/sh
# command-line to run PopLDdecay2, single thread
#
# -MaxDist 300 and the default MAF/Miss match the paper. -T 1 is the baseline the
# speed-up is measured against (PopLDdecay2 defaults to -T 4).
#
# CentOS carries two "time" programs: the one in /bin records peak memory, the
# one in /usr/bin does not. If /bin/time is missing, change /bin/time -v -p to
# time -p, but then the memory numbers are lost.

echo Start Time :
date

/bin/time	-v	-p	../../bin/PopLDdecay2	-InVCF	M5000.vcf.gz	-OutStat	M5000.stat	-MaxDist	300	-T	1
/bin/time	-v	-p	../../bin/PopLDdecay2	-InVCF	M10000.vcf.gz	-OutStat	M10000.stat	-MaxDist	300	-T	1
/bin/time	-v	-p	../../bin/PopLDdecay2	-InVCF	M20000.vcf.gz	-OutStat	M20000.stat	-MaxDist	300	-T	1
/bin/time	-v	-p	../../bin/PopLDdecay2	-InVCF	M50000.vcf.gz	-OutStat	M50000.stat	-MaxDist	300	-T	1

echo End Time :
date
