#!/bin/sh
# step4a: does PLINK use the phase that is already in the VCF?
#
# The same sites are run twice, once phased (0|1) and once unphased (0/1). A
# program that uses the phase gives different answers for the two, a program
# that re-estimates the haplotypes from the diplotype counts gives the same
# answer for both.
#
# PLINK 2.0 --r2-phased is documented to base r^2 and D' on "the maximum
# likelihood solution to the cubic equation" (Gaunt et al. 2007), i.e. on
# haplotype frequencies that are re-estimated without the phase.
#
#   sh step4a_runUnphased.sh

echo Start Time :
date

TAG=${TAG:-sub100_maf0.1}

PLINK2=${PLINK2:-$(command -v plink2 || echo /home/liaojing/01.Software/plink-260121/plink2)}
if [ -z "$PLD2" ]; then
	for c in ../../bin/PopLDdecay2 ../../../bin/PopLDdecay2; do
		if [ -x "$c" ]; then PLD2=$c; break; fi
	done
fi

IN=input/$TAG.vcf.gz
UNPH=input/${TAG}_unphased.vcf.gz
FILT="--maf 0.005 --geno 0.25 --hwe 1 1 --mind 1 --threads 8"

# bgzip and tabix come from samtools and must be in PATH:
#   source ~/01.Software/samtools-1.23/SourceMe.sh
command -v bgzip >/dev/null || { echo "ERROR: bgzip not found, source the samtools SourceMe.sh"; exit 1; }
command -v tabix >/dev/null || { echo "ERROR: tabix not found, source the samtools SourceMe.sh"; exit 1; }

mkdir -p output/plink output/poplddecay2

# --- the same sites, without the phase -------------------------------------
zcat $IN | awk 'BEGIN{OFS="\t"} /^#/{print; next} {s=$0; gsub(/\|/,"/",s); print s}' \
	| bgzip -l 9 -c > $UNPH
tabix -f -p vcf $UNPH

# --- PLINK 2.0 on the unphased copy ----------------------------------------
$PLINK2 --vcf $UNPH --make-bed --out output/plink/unph_plink2_$TAG $FILT \
	> output/plink/unph_plink2_${TAG}_makebed.log 2>&1
$PLINK2 --bfile output/plink/unph_plink2_$TAG --r2-phased cols=chrom,pos,dprimeabs,dprime \
	--ld-window-kb 300 --ld-window-r2 0 $FILT --out output/plink/unph_plink2_$TAG \
	> output/plink/unph_plink2_${TAG}_dprime.log 2>&1

# --- PopLDdecay2 on the unphased copy --------------------------------------
$PLD2 -InVCF $UNPH -OutStat output/poplddecay2/unph_$TAG -OutType 8 \
	-MaxDist 300 -T 8 > output/poplddecay2/unph_${TAG}.log 2>&1

echo "unphased runs done"

echo End Time :
date
echo "next: sh step4b_comparePhase.sh"
