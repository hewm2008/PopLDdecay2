#!/bin/sh
# make the M*.vcf.gz test datasets from the base slice
#
# Site counts 5,000 / 10,000 / 20,000 / 50,000, samples fixed at 2,504.
# Each dataset is a prefix of the one above it, so the site count is the only
# variable and the datasets stay comparable.

echo Start Time :
date

for M in 5000 10000 20000 50000; do
	zcat ../00.GetData/base_chr22_50k.vcf.gz \
		| awk -v m=$M '/^#/{print;next} {print; if(++n>=m) exit}' \
		| bgzip -c > M${M}.vcf.gz
	tabix -f -p vcf M${M}.vcf.gz
done

md5sum M*.vcf.gz > md5sum.txt
wc -l md5sum.txt

echo End Time :
date
