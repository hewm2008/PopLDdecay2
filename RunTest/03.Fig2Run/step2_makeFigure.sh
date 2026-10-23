#!/bin/sh
# make Figure 2 (needs R with ggplot2 and cowplot installed)
#
# plot_performance.r reads the .time files written by step1_runThreads.sh, prints the numbers and writes
# both the figure and perf_data.tsv into this directory.

echo Start Time :
date

Rscript	plot_performance.r

echo End Time :
date
