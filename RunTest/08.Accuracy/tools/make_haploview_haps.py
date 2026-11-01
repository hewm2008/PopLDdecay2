#!/usr/bin/env python3
"""VCF (phased) -> HapMap .haps / .legend, the input format of Haploview.

Haploview's -haps format:
  * two consecutive lines are the two chromosomes of one individual
  * every line: <familyID> <individualID> <marker1> <marker2> ...
    with one allele per marker, and the same number of columns on every line
  * the allele codes are the letters A/C/G/T, 0 for a missing call
  * the marker information file (-info) is <markername> <markerposition>,
    exactly two columns

usage: make_haploview_haps.py <in.vcf.gz> <out.haps> <out.legend>
"""
import gzip
import sys


def main():
    vcf, haps_out, legend_out = sys.argv[1], sys.argv[2], sys.argv[3]
    samples = None
    hap1, hap2 = {}, {}
    with gzip.open(vcf, "rt") as fh, open(legend_out, "w") as lg:
        for line in fh:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                hap1 = {s: [] for s in samples}
                hap2 = {s: [] for s in samples}
                continue
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            ref, alt = f[3], f[4]
            lg.write(f"{f[2]}\t{f[1]}\n")

            def base(c):
                if c == "0":
                    return ref
                if c == "1":
                    return alt
                return "0"                       # missing

            for s, g in zip(samples, f[9:]):
                gt = g.replace("|", "/")
                hap1[s].append(base(gt[0]))
                hap2[s].append(base(gt[2]))

    with open(haps_out, "w") as hp:
        for s in samples:
            hp.write(f"22\t{s}\t" + "\t".join(hap1[s]) + "\n")
            hp.write(f"22\t{s}\t" + "\t".join(hap2[s]) + "\n")
    print(f"[haps] {len(samples)} samples x 2 chromosomes -> {haps_out}")


if __name__ == "__main__":
    main()
