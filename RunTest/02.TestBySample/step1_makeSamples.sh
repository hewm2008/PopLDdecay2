#!/bin/sh
# make the N*.vcf.gz test datasets from the base slice
#
# Sample counts 100 / 250 / 500 / 2,504, sites fixed at 50,000. Each dataset is
# the first N sample columns of the one above it, so the sample count is the only
# variable. Sample order follows the original 1000G order (grouped by
# super-population), so a prefix is a full subset of all five populations.

echo Start Time :
date

for N in 100 250 500 2504; do
	zcat ../00.GetData/base_chr22_50k.vcf.gz \
		| cut -f1-9,10-$((9 + N)) \
		| bgzip -c > N${N}.vcf.gz
	tabix -f -p vcf N${N}.vcf.gz
done

md5sum N*.vcf.gz > md5sum.txt
wc -l md5sum.txt

echo End Time :
date
