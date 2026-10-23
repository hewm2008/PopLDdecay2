#!/bin/sh
# command-line to plot the LD-decay curves with the bundled perl scripts
#
# Plot_OnePop.pl draws one population per figure; Plot_MultiPop.pl overlays all
# five, which is Figure 3 in the paper. Both call Rscript internally, so R is
# needed. For the ggplot version with a vector PDF, run: Rscript plot_lddecay.r

echo Start Time :
date

/bin/time	-v	-p	perl	../../bin/Plot_OnePop.pl	-inFile	mp.AFR.stat.gz	-output	Fig3_single_AFR	-maxX	300
/bin/time	-v	-p	perl	../../bin/Plot_OnePop.pl	-inFile	mp.AMR.stat.gz	-output	Fig3_single_AMR	-maxX	300
/bin/time	-v	-p	perl	../../bin/Plot_OnePop.pl	-inFile	mp.EAS.stat.gz	-output	Fig3_single_EAS	-maxX	300
/bin/time	-v	-p	perl	../../bin/Plot_OnePop.pl	-inFile	mp.EUR.stat.gz	-output	Fig3_single_EUR	-maxX	300
/bin/time	-v	-p	perl	../../bin/Plot_OnePop.pl	-inFile	mp.SAS.stat.gz	-output	Fig3_single_SAS	-maxX	300

printf 'mp.AFR.stat.gz\tAFR\n' >  pops.list
printf 'mp.AMR.stat.gz\tAMR\n' >> pops.list
printf 'mp.EAS.stat.gz\tEAS\n' >> pops.list
printf 'mp.EUR.stat.gz\tEUR\n' >> pops.list
printf 'mp.SAS.stat.gz\tSAS\n' >> pops.list

/bin/time	-v	-p	perl	../../bin/Plot_MultiPop.pl	-inList	pops.list	-output	Fig3_subpop

echo End Time :
date
