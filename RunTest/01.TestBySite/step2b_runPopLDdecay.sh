#!/bin/sh
# command-line to run the original PopLDdecay, single thread
#
# PopLDdecay is not redistributable, so it is not in this package. step2b_runPopLDdecay.sh
# expects it next to the PopLDdecay2 directory:
#     PopLDdecay-3.45/bin/PopLDdecay   (in the package root)
# If it is missing, /bin/time reports "cannot run ... No such file or directory"
# and step3_checkOutput.sh has nothing to compare against.

echo Start Time :
date

/bin/time	-v	-p	../../PopLDdecay-3.45/bin/PopLDdecay	-InVCF	M5000.vcf.gz	-OutStat	M5000.old.stat	-MaxDist	300
/bin/time	-v	-p	../../PopLDdecay-3.45/bin/PopLDdecay	-InVCF	M10000.vcf.gz	-OutStat	M10000.old.stat	-MaxDist	300
/bin/time	-v	-p	../../PopLDdecay-3.45/bin/PopLDdecay	-InVCF	M20000.vcf.gz	-OutStat	M20000.old.stat	-MaxDist	300
/bin/time	-v	-p	../../PopLDdecay-3.45/bin/PopLDdecay	-InVCF	M50000.vcf.gz	-OutStat	M50000.old.stat	-MaxDist	300

echo End Time :
date
