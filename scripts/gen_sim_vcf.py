#!/usr/bin/env python3
"""Deterministic synthetic VCF generator for PopLDdecay2 golden-diff tests.

Produces multiple chromosomes (including a lexical-ordering trap: chr10 sorts
before chr2), variable allele-frequency waves to create LD structure, and
genotypes exercising the tricky cases (. / ., 0/., hom, het). Output is
byte-reproducible for a given --seed.
"""
import argparse
import gzip
import math
import random


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--out", required=True, help="output VCF path")
    ap.add_argument("--chroms", type=int, default=4)
    ap.add_argument("--sites", type=int, default=2000, help="sites per chromosome")
    ap.add_argument("--nsamples", type=int, default=60)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--phased", action="store_true", help="emit phased GT (a|b)")
    ap.add_argument("--gzip", action="store_true", help="gzip the output")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    bases = ["A", "C", "G", "T"]

    if args.chroms == 4:
        chroms = ["chr1", "chr2", "chr10", "chrX"]  # lexical order: chr1 < chr10 < chr2 < chrX
    else:
        chroms = [f"chr{i}" for i in range(1, args.chroms + 1)]

    samples = [f"S{i + 1}" for i in range(args.nsamples)]

    lines = [
        "##fileformat=VCFv4.2",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    ]
    for chrom in chroms:
        lines.append(f"##contig=<ID={chrom},length=1000000000>")
    lines.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + samples))

    for chrom in chroms:
        pos = 1000
        for si in range(args.sites):
            pos += rng.randint(1, 1500)
            ref = rng.choice(bases)
            alt = rng.choice([b for b in bases if b != ref])

            # allele-frequency wave -> LD; occasionally a very rare site (MAF filter)
            if si % 50 == 0:
                p = rng.uniform(0.001, 0.004)
            else:
                p = 0.3 + 0.4 * abs(math.sin(si / 37.0))

            gts = []
            for i in range(args.nsamples):
                if i == 0:
                    # first sample determines phase detection (must be a clean het)
                    gts.append("0|1" if args.phased else "0/1")
                    continue
                r = rng.random()
                if r < 0.04:
                    gts.append("./.")
                elif r < 0.08:
                    gts.append("0/.")
                else:
                    h1 = 1 if rng.random() < p else 0
                    h2 = 1 if rng.random() < p else 0
                    gts.append(f"{h1}|{h2}" if args.phased else f"{h1}/{h2}")

            lines.append("\t".join([chrom, str(pos), f"rs{si}", ref, alt, ".", "PASS", ".", "GT"] + gts))

    data = "\n".join(lines) + "\n"
    if args.gzip:
        with gzip.open(args.out, "wt") as f:
            f.write(data)
    else:
        with open(args.out, "w") as f:
            f.write(data)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()