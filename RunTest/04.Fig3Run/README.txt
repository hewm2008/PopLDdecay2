04.Fig3Run -- Figure 3: LD-decay curves of five super-populations
====================================================================

  step1a_runPopLDdecay2Multi.sh   ### step1a: one PopLDdecay2 pass over five populations
  step1b_runPopLDdecay.sh         ### step1b: the same work with the old tool, one run per population
  step1c_runPopLDdecay2Single.sh  ### step1c: the same work with the new tool, one run per population
  step2_plotCurves.sh             ### step2: draw the curves with the bundled perl scripts
  step3_checkOutput.sh            ### step3: the three designs must agree, and one pass must beat five
  plot_lddecay.r                  ### the ggplot version of the figure (vector PDF)
  sample.group                    ### sample <TAB> super-population, 2504 rows
  pop.info                        ### sample <TAB> sub-population

  cd RunTest/04.Fig3Run
  sh step1a_runPopLDdecay2Multi.sh && sh step1b_runPopLDdecay.sh \
     && sh step1c_runPopLDdecay2Single.sh \
     && sh step2_plotCurves.sh && sh step3_checkOutput.sh

Three ways to get the same five curves, which is what makes the comparison
readable:

  step1a  one indexed pass of PopLDdecay2 over the two-column sample.group:
          one .stat.gz per population (mp.AFR.stat.gz, mp.AMR.stat.gz, ...)
  step1b  the original PopLDdecay, which takes only a single-column list, so
          sample.group is split into five lists and it runs five times
  step1c  PopLDdecay2 with the same five single-column lists, five runs

step3 compares them: step1a and step1b must give the same table (the old tool
against the new one), and step1a and step1c must give the same table too (the
same tool, once against five times -- which also checks the two ways of
writing -SubPop). Then it prints how much the one-pass design saves, once
against the old tool's five runs and once against the new tool's five runs.

Note:
1. step1b needs the original PopLDdecay, which is not redistributable and not
   in this package. Put it in the package root as
   PopLDdecay-3.45/bin/PopLDdecay.

2. The input is the full chr22 (00.GetData); the 50k slice works too
   (uncomment the second input path).

3. step1a and step1c use the same -T 8, so the one-pass design is compared
   against five runs of the same tool on the same number of threads.

4. The figures need R: the perl scripts call Rscript. plot_lddecay.r is the
   ggplot version with a vector PDF; step2 uses the bundled perl scripts.
