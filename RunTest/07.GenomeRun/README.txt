07.GenomeRun -- whole-genome benchmark on 1000 Genomes Phase 3 autosomes
========================================================================

Both tools read one merged, tabix-indexed VCF of all 22 autosomes (81.3 M SNPs,
2,504 samples). This is the scale the manuscript reports on page 8; every other
directory stops at chromosome 22.

  sh step1_download.sh [local-dir]  ### 22 autosome VCFs, v5b, about 14 GB
  sh step2_merge.sh                 ### one bgzip VCF + tabix, about 14 GB, 28 min
  sh step3a_runPopLDdecay2.sh       ### -T 8, about 32 min, 21 GB peak RSS
  sh step3b_runPopLDdecay.sh        ### the original, about 61 h, 75 GB peak RSS

Measured on this machine (80 cores, 502 GB RAM, XFS):

  tool            threads   wall        peak RSS
  PopLDdecay      1         61:00:51    75 GB
  PopLDdecay2     1         3:39:37     20.6 GB
  PopLDdecay2     8         33:10       21 GB

step3b runs the original to completion, which takes about 61 hours. If you only
want to confirm that it starts, uncomment the 60-second timeout line in that
script instead: it exits 124 when the timeout fires, which is the expected
result and not a failure.

step2 is not merge-bound: `bcftools concat --naive` copies compressed blocks and
takes seconds, but `tabix` has to parse every record, and 81.3 M records of
2,504 samples is about 810 GB of VCF text once decompressed. Expect roughly
28 minutes for the index with `tabix -@ 8` (37 minutes single-threaded).

Disk: 14 GB of input plus 14 GB for the merged file; the .stat.gz outputs are a
few MB. genome/ is gitignored -- none of this ships with the package.

Note:
1. v5b is the build the paper used. The UCSC mirror carries only v5a, a
   different site set.
2. Both tools use -MaxDist 300 and -MAF 0.005, as everywhere else in RunTest.
3. step2 needs bcftools (concat --naive); bgzip and tabix come from samtools.
