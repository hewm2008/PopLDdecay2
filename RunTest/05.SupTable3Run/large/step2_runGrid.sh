#!/bin/sh
# step2: run all seven tools on every grid point of this series
#
# The steps can be run one at a time as well; this only saves typing. The bed
# built by step2c is reused by step2d and step2g, so step2c should run first.
echo Start Time :
date

for s in step2a_runPopLDdecay2.sh \
	 step2b_runPopLDdecay.sh \
	 step2c_runPLINK2.sh \
	 step2d_runPLINK19.sh \
	 step2e_runQuickLD.sh \
	 step2f_runLDkit.sh \
	 step2g_runPopLDdecay2Plink.sh; do
	echo "=== $s ==="
	sh $s || echo "WARNING: $s did not finish cleanly"
done

echo End Time :
date
