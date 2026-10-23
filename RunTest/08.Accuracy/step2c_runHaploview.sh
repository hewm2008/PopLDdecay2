#!/bin/sh
# step2c: the same pairwise LD from Haploview 4.2
#
# Haploview has no VCF reader, so the site set is converted to the HapMap
# format first (tools/make_haploview_haps.py):
#
#   input/<TAG>.vcf.gz  ->  output/haploview/<TAG>.haps + <TAG>.legend
#
# The .haps format is one line per chromosome: <familyID> <individualID> then
# one allele per marker, and two consecutive lines are the two chromosomes of
# one individual. The marker information file is <markername> <markerposition>
# and must have exactly two columns.
#
# -hwcutoff 0 excludes the markers whose HWE p-value is below 0, i.e. none of
# them: PopLDdecay2 has no HWE filter either, so all the programs see the same
# sites. The default of Haploview is 0.001, which would drop 45 of the 100 sites.
#
# Haploview.jar is downloaded on demand into input/. The Broad Institute server
# answers 403 to a plain wget, so the request carries browser headers.
#
#   sh step2c_runHaploview.sh

echo Start Time :
date

TAG=${TAG:-sub100_maf0.1}

IN=input/$TAG.vcf.gz
DIR=output/haploview
JAR=input/Haploview.jar
URL=https://www.broadinstitute.org/ftp/pub/mpg/haploview/Haploview.jar

mkdir -p $DIR

if [ ! -s $JAR ] || [ $(stat -c %s $JAR) -lt 1000000 ]; then
	echo "downloading Haploview.jar"
	curl -sSL -o $JAR \
		-H 'User-Agent: Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/125.0.0.0 Safari/537.36' \
		-H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,*/*;q=0.8' \
		-H 'Accept-Language: en-US,en;q=0.9' \
		-H 'Accept-Encoding: gzip, deflate, br' \
		-H 'Referer: https://www.broadinstitute.org/software/haploview/haploview' \
		-H 'Sec-Fetch-Dest: document' -H 'Sec-Fetch-Mode: navigate' -H 'Sec-Fetch-Site: cross-site' \
		$URL
fi
echo "Haploview.jar: $(stat -c %s $JAR) bytes"
if [ $(stat -c %s $JAR) -lt 1000000 ]; then
	echo "ERROR: Haploview.jar was not downloaded (the server answers 403 here)."
	echo "       Fetch it with a browser and put it into $JAR, then run this step again."
	exit 1
fi

python3 tools/make_haploview_haps.py $IN $DIR/$TAG.haps $DIR/$TAG.legend

java -Xmx32g -jar $JAR -nogui -memory 16384 \
	-dprime -out $DIR/$TAG \
	-haps $DIR/$TAG.haps -info $DIR/$TAG.legend \
	-maxdistance 300 -minMAF 0.005 -hwcutoff 0 -missingCutoff 1.0 \
	> $DIR/$TAG.log 2>&1

echo "haploview pairs: $(( $(wc -l < $DIR/$TAG.LD) - 1 ))"

echo End Time :
date
echo "next: sh step3_compare.sh"
