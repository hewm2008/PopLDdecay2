#!/bin/sh
# step1: the first N sites with MAF >= MAF
#
#   ../00.GetData/base_chr22_50k.vcf.gz  ->  input/sub<N>_maf<M>.vcf.gz
#
# Why not simply the first N sites of the slice: at MAF 0.005 nearly every R^2
# is 0.000xx and the comparison says nothing. At MAF 0.1 the sites are common
# and the LD between them is strong enough to be informative.
#
# The sites are picked the same way 00.GetData/step4_makeSubset.sh does it:
# PopLDdecay2 selects the sites that pass -MAF / -Miss, and they are then taken
# back out of the VCF with their genotypes.
#
# The site count and the MAF are written into this script; change them here and
# the name of the file follows, so several sets can live side by side:
#
#   NSITE=100   MAF=0.1   ->  input/sub100_maf0.1.vcf.gz

echo Start Time :
date

SRC=../00.GetData/base_chr22_50k.vcf.gz
NSITE=100
MAF=0.1
TAG=sub${NSITE}_maf${MAF}

# the binary sits in bin/ of the package root; override with PLD2=... if not
if [ -z "$PLD2" ]; then
	for c in ../../bin/PopLDdecay2 ../../../bin/PopLDdecay2; do
		if [ -x "$c" ]; then PLD2=$c; break; fi
	done
fi
echo "PopLDdecay2: $PLD2"

# bgzip and tabix come from samtools and must be in PATH:
#   source ~/01.Software/samtools-1.23/SourceMe.sh
command -v bgzip >/dev/null || { echo "ERROR: bgzip not found, source the samtools SourceMe.sh"; exit 1; }
command -v tabix >/dev/null || { echo "ERROR: tabix not found, source the samtools SourceMe.sh"; exit 1; }

mkdir -p input

W=$(mktemp -d)
trap 'rm -rf $W' EXIT

echo "1/3 sites that pass -MAF $MAF / -Miss 0.25"
$PLD2 -InVCF $SRC -OutStat $W/ref -MaxDist 300 -T 16 -MAF $MAF -OutFilterSNP >/dev/null 2>&1
zcat $W/ref.vcf.filter.gz | grep -v '^#' | cut -f2 > $W/pass.pos
echo "    $(wc -l < $W/pass.pos) sites pass"

echo "2/3 take the first $NSITE of them"
head -$NSITE $W/pass.pos > $W/want.pos

echo "3/3 pull them back out of the VCF, bi-allelic SNPs only"
zcat $SRC | awk -v m=$NSITE '
	BEGIN{OFS="\t"}
	NR==FNR { want[$1]=1; next }
	/^#/ { print; next }
	($2 in want) && length($4)==1 && length($5)==1 &&
	$4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/ && $4 != $5 {
		$3=$2; print; if (++n >= m) exit
	}' $W/want.pos - \
	| bgzip -l 9 -c > input/$TAG.vcf.gz
tabix -f -p vcf input/$TAG.vcf.gz

echo "sites  : $(zcat input/$TAG.vcf.gz | grep -vc '^#')"
echo "samples: $(zcat input/$TAG.vcf.gz | grep -m1 '^#CHROM' | awk '{print NF-9}')"
echo "range  : $(zcat input/$TAG.vcf.gz | grep -v '^#' | cut -f2 | head -1) - $(zcat input/$TAG.vcf.gz | grep -v '^#' | cut -f2 | tail -1)"

echo End Time :
date
echo "next: sh step2a_runPopLDdecay2.sh"
