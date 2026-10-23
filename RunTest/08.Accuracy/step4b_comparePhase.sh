#!/bin/sh
# step4b: report the phase test
#
#   tools/phase_test.py   compares the phased and the unphased run of each
#                         program and says whether the phase was used
#
#   sh step4b_comparePhase.sh

echo Start Time :
date

TAG=${TAG:-sub100_maf0.1}

python3 tools/phase_test.py $TAG

echo End Time :
date
echo "done. read output/compare/$TAG/report.txt and .../phase.txt"
