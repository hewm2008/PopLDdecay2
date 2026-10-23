The latest version of PopLDdecay2, scripts and datasets generated in this manuscript were uploaded.

- 00.GetData         ## prepare the data: download chr22, rebuild the index, cut the slice
- 01.TestBySite      ## benchmark with the increase of the number of SNP sites
- 02.TestBySample    ## benchmark with the increase of the number of samples
- 03.Fig2Run         ## Figure 2: real time, peak RSS, strong-scaling efficiency
- 04.Fig3Run         ## Figure 3: LD-decay curves of five super-populations
                      ## (one pass, five runs of the old tool, five runs of the new one)
- 05.SupTable3Run    ## Supplementary Table S3: matched benchmark against
                      ## PopLDdecay 3.45 / PLINK 2 / PLINK 1.9 / quickLD / LDkit,
- 06.Verify          ## the regression gate for the binary you just built
- 07.GenomeRun       ## whole-genome run: all 22 autosomes, 81.3 M SNPs
- 08.Accuracy        ## the pairwise LD of PopLDdecay2 against PLINK 2.0 and
                      ## Haploview, on 100 sites with the true haplotypes as the
                      ## ground truth

Note:
1. The .tbi index that comes with the chr22 VCF from 1000 Genomes is STALE: the
   data was rewritten after the index was built. htslib prints
   "[W::hts_idx_load3] The index file is older than the data file", which is easy
   to miss. With a stale index the multi-threaded fast path silently turns off,
   so -T 8 runs as slowly as -T 1. Always run 00.GetData/step3_makeTBI.sh before
   benchmarking.

2. PopLDdecay (the original) is not redistributable and is not included here.
   The step2b / step3b scripts expect it in the package root:
       PopLDdecay-3.45/bin/PopLDdecay
   If it is absent, /bin/time reports "cannot run" and the check steps report
   the missing .old.stat.gz files.

3. The tests were run on CentOS Stream, which carries two "time" programs: the
   one in /bin records detailed information including peak memory, while the
   one in /usr/bin records only the running time. Please check that /bin/time
   exists; if it does not, change "/bin/time -v -p" to "time -p" -- but then the
   peak-memory numbers of Figure 2 cannot be reproduced.

4. Data: 1000 Genomes Project Phase 3 (20130502) chr22, 2,504 samples.
   The EBI copy is v5b; the UCSC mirror only carries v5a, a different site set.
