#!/bin/sh
# step2b: the same pairwise LD from PLINK 2.0
#
#   plink2 --r2-phased cols=chrom,pos,dprimeabs,dprime
#
# PLINK 1.9 --r2 dprime was run on the same data and returned bit-identical
# numbers (Pearson 1.000000, max|d| 0 over all site pairs), so only PLINK 2.0
# is run here.
#
# The filters are aligned with PopLDdecay2: -MAF 0.005, -Miss 0.25, and NO HWE
# filter, because PopLDdecay2 has none either. PLINK 2.0 excludes variants with
# a HWE p-value below p * 10^(-n*k), so "--hwe 1 1" gives a threshold of
# 10^-2504 and removes nothing (with --hwe 1e-30 it would drop 45 of the 100
# sites, and the comparison would shrink to a quarter of the site pairs).
#
#   sh step2b_runPLINK2.sh

echo Start Time :
date

TAG=${TAG:-sub100_maf0.1}

PLINK2=${PLINK2:-$(command -v plink2 || echo /home/liaojing/01.Software/plink-260121/plink2)}
echo "PLINK 2.0: $PLINK2"

IN=input/$TAG.vcf.gz
DIR=output/plink

mkdir -p $DIR

FILT="--maf 0.005 --geno 0.25 --hwe 1 1 --mind 1 --threads 8"

$PLINK2 --vcf $IN --make-bed --out $DIR/plink2_$TAG $FILT > $DIR/plink2_${TAG}_makebed.log 2>&1
$PLINK2 --bfile $DIR/plink2_$TAG --r2-phased cols=chrom,pos,dprimeabs,dprime \
	--ld-window-kb 300 --ld-window-r2 0 $FILT --out $DIR/plink2_$TAG \
	> $DIR/plink2_${TAG}_dprime.log 2>&1

echo "plink2 sites : $(grep -c . $DIR/plink2_$TAG.bim)"
echo "plink2 pairs : $(( $(wc -l < $DIR/plink2_$TAG.vcor) - 1 ))"

echo End Time :
date
echo "next: sh step3_compare.sh"
