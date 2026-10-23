#!/bin/sh
# step2a: run PopLDdecay2 with -OutType 8 on the site set
#
# -OutType 8 writes the pairwise LD table (D' / LOD / r^2 / CIlow / CIhi / Dist)
# to <prefix>.LD.gz and the distance-binned LD decay statistics to <prefix>.stat.gz.
# The other programs are compared against the .LD.gz table in step3.
#
#   sh step2a_runPopLDdecay2.sh

echo Start Time :
date

TAG=${TAG:-sub100_maf0.1}

if [ -z "$PLD2" ]; then
	for c in ../../bin/PopLDdecay2 ../../../bin/PopLDdecay2; do
		if [ -x "$c" ]; then PLD2=$c; break; fi
	done
fi
echo "PopLDdecay2: $PLD2"

IN=input/$TAG.vcf.gz
OUT=output/poplddecay2/$TAG

mkdir -p output/poplddecay2

$PLD2 -InVCF $IN -OutStat $OUT -OutType 8 -MaxDist 300 -T 8 -OutFilterSNP

echo "pairs: $(( $(zcat $OUT.LD.gz | wc -l) - 1 ))"
echo "sites kept by -MAF/-Het/-Miss: $(zcat $OUT.vcf.filter.gz | grep -vc '^#')"

echo End Time :
date
echo "next: sh step3_compare.sh"
