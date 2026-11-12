#!/bin/sh
# step2: merge the 22 autosomes into one bgzip VCF and index it (about 14 GB)
#
# bcftools concat --naive appends the files in order without recompressing and
# checks that the headers agree; no input index is needed. The output is bgzip,
# so tabix can index it and PopLDdecay2's fast path can use it.
#
# The order is numeric, not the lexicographic order of a glob (chr1, chr10,
# chr11, ..., chr2), which concat would reject as unsorted.
#
# Indexing is the slow part, not the concat: tabix has to parse every record,
# and 81 M records of 2,504 samples is about 810 GB of VCF text. -@ 8 buys
# roughly a quarter of the time back; the cost is parse-bound, not I/O-bound.

P=phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
G=genome/ALL.autosomes.phase3.vcf.gz

FILES=""
for i in $(seq 1 22); do
	FILES="$FILES genome/ALL.chr${i}.$P"
done

echo Start Time :
date

bcftools concat --naive -O z -o $G $FILES
tabix -@ 8 -f -p vcf $G

ls -lh $G
tabix -l $G | tr '\n' ' '

echo End Time :
date
