#!/usr/bin/env python3
"""Generate a PopLDdecay `-InGenotype` input file from a biallelic VCF.

Format: `chr<TAB>site<TAB>IUPAC1 IUPAC2 ... IUPACn` (one IUPAC code per
sample; A/C/G/T = homozygote, M/K/Y/R/W/S = heterozygote, '-' = missing),
mirroring the VCF2Genotype -WithHeader -NoRef style the original tool expects.

An optional `#CHROM  POS <S1> <S2> ...` header is written when --with-header
is given (needed for `-SubPop` runs).

IUPAC het code chosen from the site's REF/ALT bases; sample GT 0/1|1/0 -> het,
0/0 -> REF homo, 1/1 -> ALT homo, './.' -> '-'.
"""
import argparse
import gzip
import random
import sys

HET = {
    frozenset("AC"): "M", frozenset("GT"): "K", frozenset("CT"): "Y",
    frozenset("AG"): "R", frozenset("AT"): "W", frozenset("CG"): "S",
}


def iupac_pair(ref, alt):
    for bases, code in HET.items():
        if bases == frozenset((ref, alt)):
            return code
    return "N"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--in-vcf", required=True)
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--with-header", action="store_true")
    ap.add_argument("--nsamples", type=int, default=0)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    openf = gzip.open if args.in_vcf.endswith(".gz") else open
    all_samples = []
    with openf(args.in_vcf, "rt") as f:
        samples = []
        lines = []
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                toks = line.rstrip("\n").split("\t")
                all_samples = toks[9:]
                samples = all_samples[:]
                if args.nsamples:
                    rng = random.Random(args.seed)
                    samples = rng.sample(all_samples, min(args.nsamples, len(all_samples)))
                continue
            lines.append(line.rstrip("\n"))

    out = open(args.out, "w")
    if args.with_header:
        out.write("#CHROM\tPOS\t" + " ".join(samples) + "\n")

    for line in lines:
        toks = line.split("\t")
        chrom, pos, ref, alt = toks[0], toks[1], toks[3].upper(), toks[4].upper()
        if len(ref) != 1 or len(alt) != 1:
            continue
        gts = toks[9:]
        gtmap = []
        for s in samples:
            idx = all_samples.index(s)
            g = gts[idx]
            a1 = g.split(":")[0]
            if "/" not in a1 and "|" not in a1:
                gtmap.append("-")
                continue
            sep = "/" if "/" in a1 else "|"
            x, y = a1.split(sep)
            if x == "." and y == ".":
                gtmap.append("-")
            elif x == "0" and y == "0":
                gtmap.append(ref)
            elif x == "1" and y == "1":
                gtmap.append(alt)
            elif (x, y) in (("0", "1"), ("1", "0")):
                gtmap.append(iupac_pair(ref, alt))
            else:
                gtmap.append("-")
        out.write("%s\t%s\t%s\n" % (chrom, pos, " ".join(gtmap)))
    out.close()
    print("wrote %d sites, %d samples" % (len(lines), len(samples)))


if __name__ == "__main__":
    sys.exit(main())
