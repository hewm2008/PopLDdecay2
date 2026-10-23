00.GetData -- prepare the data
================================

  step1_download.sh    ### download the 1000 Genomes Phase 3 chr22 VCF and the sample panel
  step2_makeGroup.sh   ### make sample.group (sample <TAB> super-population) from the panel
  step3_makeTBI.sh     ### rebuild the tabix index -- DO NOT SKIP, the downloaded one is stale
  step4_makeSubset.sh  ### make the 50k slice that ships with this package
  base_chr22_50k.vcf.gz  ### 50,000 bi-allelic SNPs x 2,504 samples (19 MB)
  base_chr22_50k.vcf.gz.tbi  ### its index, ready to use
  chr22/               ### the full chr22, once downloaded (206 MB)

In order:
  sh step1_download.sh        # only if you do not have the data already
  sh step2_makeGroup.sh
  sh step3_makeTBI.sh
  sh step4_makeSubset.sh      # normally not needed: the slice ships with the package

Note:
1. The index that comes with the chr22 download is STALE: the data was rewritten
   after the index was built. htslib prints
   "[W::hts_idx_load3] The index file is older than the data file", which is easy
   to miss. With a stale index the multi-threaded fast path silently turns off,
   so -T 8 runs as slowly as -T 1 and the speed-up is badly underestimated.
   Always run step3_makeTBI.sh before benchmarking.

2. The 50k slice is enough for 01-04 and 08: those need no download at all.
   The full chr22 is needed for 05 and 07; 03 and 04 use it by default but also
   run on the 50k slice (uncomment the second input path).

3. The slice is cut from the sites that pass PopLDdecay2's own filters, so it
   holds nothing but bi-allelic SNPs and its site count is exact (see the
   comment in step4_makeSubset.sh).

4. Data: 1000 Genomes Project Phase 3 (20130502) chr22, 2,504 samples. The EBI
   copy is v5b; the UCSC mirror only carries v5a, a different site set.
