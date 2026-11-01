01.TestBySite -- how the runtime scales with the number of sites
=============================================================

  step1_makeSites.sh          ### step1: make M5000 / M10000 / M20000 / M50000 from the slice
  step2a_runPopLDdecay2.sh    ### step2a: run PopLDdecay2, -T 1, timed
  step2b_runPopLDdecay.sh     ### step2b: run the original PopLDdecay, the baseline
  step3_checkOutput.sh        ### step3: the two outputs must be byte-identical (diff = 0)
  md5sum.txt                  ### the digests of the generated datasets

  cd RunTest/01.TestBySite
  sh step1_makeSites.sh && sh step2a_runPopLDdecay2.sh \
                       && sh step2b_runPopLDdecay.sh && sh step3_checkOutput.sh

The site counts are 5,000 / 10,000 / 20,000 / 50,000 and the samples are fixed
at 2,504. Each dataset is a prefix of the next one, so the site count is the
only variable and the runs stay comparable.

Note:
1. step2b needs the original PopLDdecay, which is not redistributable and not
   in this package. Put it in the package root as
   PopLDdecay-3.45/bin/PopLDdecay. Without it /bin/time reports "cannot run"
   and step3 has nothing to compare against.

2. The .stat.gz containers differ (BGZF vs plain gzip), so the file sizes
   differ; only the decompressed content is compared.

3. /bin/time -v records the peak memory. If /bin/time is missing, change it to
   time -p -- but then the memory numbers are lost.
