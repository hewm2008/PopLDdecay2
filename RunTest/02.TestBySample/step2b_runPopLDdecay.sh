#!/bin/sh
# command-line to run the original PopLDdecay, single thread
#
# PopLDdecay is not redistributable, so it is not in this package. step2b_runPopLDdecay.sh
# expects it next to the PopLDdecay2 directory:
#     PopLDdecay-3.45/bin/PopLDdecay   (in the package root)

echo Start Time :
date

/bin/time	-v	-p	../../PopLDdecay-3.45/bin/PopLDdecay	-InVCF	N100.vcf.gz	-OutStat	N100.old.stat	-MaxDist	300
/bin/time	-v	-p	../../PopLDdecay-3.45/bin/PopLDdecay	-InVCF	N250.vcf.gz	-OutStat	N250.old.stat	-MaxDist	300
/bin/time	-v	-p	../../PopLDdecay-3.45/bin/PopLDdecay	-InVCF	N500.vcf.gz	-OutStat	N500.old.stat	-MaxDist	300
/bin/time	-v	-p	../../PopLDdecay-3.45/bin/PopLDdecay	-InVCF	N2504.vcf.gz	-OutStat	N2504.old.stat	-MaxDist	300

echo End Time :
date
