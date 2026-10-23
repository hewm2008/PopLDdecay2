08.Accuracy -- the accuracy of the pairwise LD of PopLDdecay2, against PLINK 2.0 and Haploview
===============================================================

  step1_makeSubset.sh        ### step1: the first 100 sites with MAF >= 0.1
  step2a_runPopLDdecay2.sh   ### step2a: PopLDdecay2 -OutType 8
  step2b_runPLINK2.sh        ### step2b: PLINK 2.0 --r2-phased
  step2c_runHaploview.sh     ### step2c: Haploview 4.2 -nogui -dprime
  step3a_compare.sh          ### step3a: the R^2 of every program against the ground truth
  step3b_showPairs.sh        ### step3b: the three tables side by side, a quick look
  step4a_runUnphased.sh      ### step4a: the same sites without the phase, run again
  step4b_comparePhase.sh     ### step4b: did the program use the phase?

  cd RunTest/08.Accuracy
  sh step1_makeSubset.sh && sh step2a_runPopLDdecay2.sh && sh step2b_runPLINK2.sh \
     && sh step2c_runHaploview.sh && sh step3a_compare.sh && sh step3b_showPairs.sh \
     && sh step4a_runUnphased.sh && sh step4b_comparePhase.sh

step2a / step2b / step2c only read the same input, so they can run at the same
time:  sh step2a_runPopLDdecay2.sh & sh step2b_runPLINK2.sh & sh step2c_runHaploview.sh & wait

  input/   sub100_maf0.1.vcf.gz, Haploview.jar
  output/  poplddecay2/ plink/ haploview/ and compare/sub100_maf0.1/ (report.txt,
           the tables, phase.txt)
  tools/   ld_truth.py compare_ld.py phase_test.py make_haploview_haps.py

The 100 sites are the first ones of ../00.GetData/base_chr22_50k.vcf.gz with
MAF >= 0.1 (2,504 samples, phased); they give 3,590 site pairs inside the
300 kb window. The ground truth is the phased haplotypes of the VCF, counted
directly by tools/ld_truth.py. step3a compares the three programs against it;
step3b only puts their tables next to each other (the pairs come in the same
order, so the rows line up).

PLINK does not use the phase that is already in the VCF: step4 runs the same
sites phased and unphased and PLINK returns identical numbers for both, while
PopLDdecay2 and Haploview use the phase and agree with the truth.

Note:
1. bgzip and tabix come from samtools and must be in PATH:
       source ~/01.Software/samtools-1.23/SourceMe.sh
2. Haploview has no VCF reader; tools/make_haploview_haps.py writes the HapMap
   format it wants. step2c downloads Haploview.jar into input/ (the server
   answers 403 to a plain wget, so the request carries browser headers).
3. The scripts use paths relative to this directory, so run them from here; the
   tag can be changed with TAG=<name>. The screen gets the headline table and a
   short note (step3a) and the two verdicts (step4b); report.txt and phase.txt
   have everything (--full prints it).
