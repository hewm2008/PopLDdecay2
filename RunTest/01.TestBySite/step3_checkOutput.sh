#!/bin/sh
# the two outputs must be byte-identical after gunzip (diff = 0)
#
# The .stat.gz containers differ (BGZF vs plain gzip), so the file sizes differ;
# only the decompressed content is compared. Needs step2a_runPopLDdecay2.sh and step2b_runPopLDdecay.sh first.

for f in M5000 M10000 M20000 M50000; do
	[ -f $f.stat.gz ] || { echo "$f  no $f.stat.gz -- run step2a_runPopLDdecay2.sh first"; continue; }
	[ -f $f.old.stat.gz ] || { echo "$f  no $f.old.stat.gz -- run step2b_runPopLDdecay.sh first"; continue; }
	zcat $f.old.stat.gz > out.old
	zcat $f.stat.gz     > out.new
	echo "$f	diff=$(diff out.old out.new | wc -l)"
done
rm -f out.old out.new
