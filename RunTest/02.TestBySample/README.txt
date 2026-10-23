02.TestBySample -- how the runtime scales with the number of samples
===================================================================

  step1_makeSamples.sh        ### step1: make N100 / N250 / N500 / N2504 from the slice
  step2a_runPopLDdecay2.sh    ### step2a: run PopLDdecay2, -T 1, timed
  step2b_runPopLDdecay.sh     ### step2b: run the original PopLDdecay, the baseline
  step3_checkOutput.sh        ### step3: the two outputs must be byte-identical (diff = 0)
  md5sum.txt                  ### the digests of the generated datasets

  cd RunTest/02.TestBySample
  sh step1_makeSamples.sh && sh step2a_runPopLDdecay2.sh \
                        && sh step2b_runPopLDdecay.sh && sh step3_checkOutput.sh

The sample counts are 100 / 250 / 500 / 2,504 and the sites are fixed at 50,000.
Each dataset is the first N sample columns of the next one, so the sample count
is the only variable. The sample order follows the original 1000G order, so a
prefix is a full subset of all five super-populations.

Note:
1. step2b needs the original PopLDdecay, which is not redistributable and not
   in this package. Put it in the package root as
   PopLDdecay-3.45/bin/PopLDdecay.

2. The .stat.gz containers differ (BGZF vs plain gzip), so the file sizes
   differ; only the decompressed content is compared.

3. This is where the memory scaling shows: the sites are held as 2-bit
   genotypes, so the cost per site grows with the sample count.
