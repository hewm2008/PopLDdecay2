Supplementary Table S3: matched benchmark against PopLDdecay 3.45, PLINK 2,
PLINK 1.9, quickLD and LDkit
=====================================================================

Two series, same four sample counts (250 / 500 / 1000 / 2504):

  small/   ### sites 2,000 / 5,000 / 10,000 / 50,000              (16 cells;
           ### 100,000 is used only from large/ to avoid a duplicate)
  large/   ### sites 100,000 / 200,000 / 300,000 / 500,000       (16 cells,
           ### the grid of the published table)

  cd small   &&  sh step1_makeGrid.sh  &&  sh step2_runGrid.sh
  cd large   &&  sh step1_makeGrid.sh  &&  sh step2_runGrid.sh

  sh run_all.sh  ### V3: 32 cells x 5 repeats, quickLD last

Every step accepts SITES= and SAMPLES= overrides, e.g. a single cell:
  SITES=2000 SAMPLES=250 sh step1_makeGrid.sh && SITES=2000 SAMPLES=250 sh step2_runGrid.sh

  step1_makeGrid.sh              ### step1: make input/ for this series
  step2a_runPopLDdecay2.sh       ### PopLDdecay2 on the VCF, -T 1 and -T 8
  step2b_runPopLDdecay.sh        ### PopLDdecay 3.45, the baseline (no -T)
  step2c_runPLINK2.sh            ### PLINK 2, VCF -> bed, --r2-phased, --threads 8
  step2d_runPLINK19.sh           ### PLINK 1.9, --r2 square, --threads 8 (memory-guarded)
  step2e_runQuickLD.sh           ### quickLD, 8 threads, timeout, memory-guarded
  step2f_runLDkit.sh             ### LDkit, --threads 8
  step2g_runPopLDdecay2Plink.sh  ### PopLDdecay2 -InPlink on the same bed as PLINK
  step2_runGrid.sh               ### runs step2a..step2g in this order
  step3_makeSummary.py           ### mean of 5 repeats -> s3_summaryV3.tsv

Every tool reads the SAME input/n<N>_m<M>.vcf.gz (or the one bed made from it
by step2c), with MAF 0.005, a 300 kb window and r^2, so the timings are
comparable. Outputs go to out/<tool>/ next to .r1..r5.time files holding the
/bin/time -v records; s3_summaryV3.tsv reports the mean of the 5 repeats.
run_all.sh uses 2,000/5,000/10,000/50,000 in small/ and
100,000/200,000/300,000/500,000 in large/; 100k is kept only once (large).

The grid. Each cell is a systematic sample of the biallelic SNPs of chr22:
record i of the T = 1,055,454 biallelic SNPs is kept when (i*M) mod T < M,
which is true exactly M times, so the sites are spread over the whole
chromosome (not taken from its start) and every position shard of the -T 8
run has work. Only biallelic SNPs are used: PLINK 2 refuses a .bim with
multiallelic variants and the other tools skip them, so a mixed grid would
not give the tools the same sites.

Tool paths. Each script names its tool at the top and every path can be
overridden from the environment:
  PopLDdecay2     ../../../bin/PopLDdecay2
  PopLDdecay 3.45 ../../../PopLDdecay-3.45/bin/PopLDdecay   (not shipped, see
                  RunTest/README.txt note 2)
  PLINK 2/1.9     from PATH, else the copies used for the paper
  LDkit, quickLD  ../../../rebuttal/tools/...  (not part of this package:
                  set LDKIT_JAR / QUICKLD_DIR, or copy them there)

Note:
1. quickLD is run with -threads 8 up to 200k sites, plus a timeout
   (QLD_TIMEOUT, default 4 h); 300k/500k are skipped as too large. In this
   release -threads computes no pairs, so those cells are marked EMPTY_OUTPUT
   (or SKIPPED_MEMORY / SKIPPED_SIZE) in s3_summaryV3.tsv. quickLD also has no
   distance window on the correlation -- it emits every pair of the parsed
   variants -- so its column is an upper bound and not the same work as the
   others.

2. PLINK 1.9 --r2 square holds a kept x kept matrix: memory and output grow
   with the square of the kept variant count (about 3.3 GB at 20,700 kept
   variants, 80 GB at 102,600). Cells whose projection does not fit in three
   times the free memory are skipped with the projection printed. The .ld.gz
   is deleted after each repeat; only the .time records are kept.

3. The VCF -> bed conversion of the PLINKs (and of -InPlink) is not part of
   the timed LD step; its record is kept in out/bed/<cell>.convert.time. The
   PLINK columns are therefore a lower bound of the end-to-end cost.

4. LDkit deletes the tabix index of the file it reads, and that index is what
   makes PopLDdecay2 -T 8 fast. step2f rebuilds it after every run; step2a and
   step2e rebuild it too if it is missing.

5. -InPlink has no indexed fast path and no position sharding, so -T there
   only parallelises the decoding. The bed it reads is the one PLINK 2 reads,
   which is what makes the two engines comparable without a format cost. A
   .bed stores unphased genotypes, so this step computes the unphased r^2 --
   unlike the phased VCF runs of step2a.

6. quickLD and PLINK write one line per SNP pair; in the large series this is
   where the disk goes. Delete out/ after collecting the .time records.
