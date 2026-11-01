03.Fig2Run -- Figure 2: real time, peak RSS, strong-scaling efficiency
=====================================================================

  step1_runThreads.sh   ### step1: run PopLDdecay2 with -T 1/4/8/16/32, five repeats each
  step2_makeFigure.sh   ### step2: Rscript plot_performance.r -- reads the .time files
  plot_performance.r    ### the figure script (R with ggplot2 and cowplot)
  perf_data.tsv         ### the numbers behind the figure, written by plot_performance.r

  cd RunTest/03.Fig2Run
  sh step1_runThreads.sh && sh step2_makeFigure.sh

The figure has three panels: run time, peak RSS and strong-scaling efficiency
(E = T1_wall / (T * wall_at_T)). The five repeats are what the error bars are
made of -- do not cut them to one. The original PopLDdecay has no -T and is run
once per repeat as the baseline.

Note:
1. step1 needs the full chr22 (00.GetData). Uncomment the second input path to
   run on the 50k slice instead: the absolute times are much smaller, the shape
   of the curves is the same.

2. Every run writes its /bin/time -v record to a .time file; plot_performance.r
   reads those, so nothing is typed into the figure by hand. If /bin/time is
   missing, change /bin/time -v -p to time -p -- but then the memory panel
   cannot be reproduced.

3. The figure needs R with ggplot2 and cowplot installed.
