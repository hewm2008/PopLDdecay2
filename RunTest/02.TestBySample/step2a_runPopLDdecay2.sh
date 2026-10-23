#!/bin/sh
# command-line to run PopLDdecay2, single thread
#
# Sites fixed at 50,000; samples 100 / 250 / 500 / 2,504. -MaxDist 300 and the
# default MAF/Miss match the paper. CentOS carries two "time" programs: the one
# in /bin records peak memory, the one in /usr/bin does not.

echo Start Time :
date

/bin/time	-v	-p	../../bin/PopLDdecay2	-InVCF	N100.vcf.gz	-OutStat	N100.stat	-MaxDist	300	-T	1
/bin/time	-v	-p	../../bin/PopLDdecay2	-InVCF	N250.vcf.gz	-OutStat	N250.stat	-MaxDist	300	-T	1
/bin/time	-v	-p	../../bin/PopLDdecay2	-InVCF	N500.vcf.gz	-OutStat	N500.stat	-MaxDist	300	-T	1
/bin/time	-v	-p	../../bin/PopLDdecay2	-InVCF	N2504.vcf.gz	-OutStat	N2504.stat	-MaxDist	300	-T	1

echo End Time :
date
