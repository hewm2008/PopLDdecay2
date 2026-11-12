#!/bin/sh
# the three designs must give identical curves, and one pass must beat five
#
# diff compares the decompressed content (the containers differ: BGZF vs gzip).
#   step1a  one PopLDdecay2 pass over five populations (two-column -SubPop)
#   step1b  the old tool, five runs (one run per population)
#   step1c  PopLDdecay2, five runs (one run per population)
# The speed-up of the one-pass design is measured against both.

for p in AFR AMR EAS EUR SAS; do
	[ -f mp.$p.stat.gz ] || { echo "$p  no mp.$p.stat.gz -- run step1a_runPopLDdecay2Multi.sh first"; continue; }
	[ -f $p.old.stat.gz ] || { echo "$p  no $p.old.stat.gz -- run step1b_runPopLDdecay.sh first"; continue; }
	zcat $p.old.stat.gz > out.old
	zcat mp.$p.stat.gz  > out.new
	echo "$p	diff=$(diff out.old out.new | wc -l)"
done
rm -f out.old out.new

# The same tool, five runs instead of one: the sites and the pairs are the same,
# so the tables must be byte-identical. This also checks the two ways of
# writing -SubPop (one two-column file against five single-column files).
for p in AFR AMR EAS EUR SAS; do
	[ -f $p.NewS.stat.gz ] || { echo "$p  no $p.NewS.stat.gz -- run step1c_runPopLDdecay2Single.sh first"; continue; }
	zcat $p.NewS.stat.gz > out.single
	zcat mp.$p.stat.gz   > out.pass
	echo "$p	single-vs-pass	diff=$(diff out.single out.pass | wc -l)"
done
rm -f out.single out.pass

# "Elapsed (wall clock) time (h:mm:ss or m:ss): 1:23.45" -> seconds
el() { grep Elapsed $1 | awk '{print $NF}' | awk -F: '{print (NF==3)? $1*3600+$2*60+$3 : $1*60+$2}'; }

[ -f mp.time ] || exit 0
new=$(el mp.time)
old=0
have_old=1
for p in AFR AMR EAS EUR SAS; do
	[ -f $p.old.time ] || { have_old=0; break; }
	old=$(awk -v a=$old -v b=$(el $p.old.time) 'BEGIN{print a+b}')
done
if [ "$have_old" = 1 ]; then
	awk -v o=$old -v n=$new 'BEGIN{printf "one pass %.1fs vs five runs of the old tool %.1fs = %.1fx\n", n, o, o/n}'
fi

# The same comparison with the new tool on both sides: what the one-pass
# design itself saves, apart from the tool being faster.
single=0
have_single=1
for p in AFR AMR EAS EUR SAS; do
	[ -f $p.NewS.time ] || { have_single=0; break; }
	single=$(awk -v a=$single -v b=$(el $p.NewS.time) 'BEGIN{print a+b}')
done
if [ "$have_single" = 1 ]; then
	awk -v o=$single -v n=$new 'BEGIN{printf "one pass %.1fs vs five runs of the new tool %.1fs = %.1fx\n", n, o, o/n}'
fi
