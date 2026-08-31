#!/bin/bash
echo Start Time : 
date
##   download the really data  ###
#wget  -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz 
#wget  -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
# cut -f  1,3   integrated_call_samples_v3.20130502.ALL.panel   > sample.group ; gzip  sample.group
# re-tabix ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz by  tabix # 
echo End Time : 
date
