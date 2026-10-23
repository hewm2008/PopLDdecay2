#!/bin/sh
# step3: rebuild the tabix index -- DO NOT SKIP
#
# The .tbi that comes with the chr22 VCF from 1000 Genomes is stale: the data was
# rewritten after the index was built. htslib prints
#   [W::hts_idx_load3] The index file is older than the data file
# and this is easy to miss. With a stale index the multi-threaded fast path
# silently turns off, so -T 8 runs as slowly as -T 1 and the speed-up is badly
# underestimated. Rebuilding takes about 30 seconds.

echo Start Time :
date

rm -f chr22/chr22.vcf.gz.tbi
tabix -f -p vcf chr22/chr22.vcf.gz

# should list records, with no "index file is older" warning
tabix chr22/chr22.vcf.gz 22:16050000-16100000 | head -2 | cut -f1-5

echo End Time :
date
