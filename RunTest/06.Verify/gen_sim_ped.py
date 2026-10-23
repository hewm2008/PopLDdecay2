#!/usr/bin/env python3
"""Generate a PED/MAP pair for `bin/mis/plink2genotype.pl`, plus the `.geno`
that the converter is expected to produce from it.

The three files come from one allele matrix, exactly as
`scripts/gen_sim_plink.py` does for the `.bed` path: the converter's output is
compared byte-for-byte against the `.geno` written here, so a change in the
converter cannot pass unnoticed. Until now this script had no test at all, in
either the golden suite or the CI gate.

Outputs, for `-o PREFIX`:
  PREFIX.ped    classic PED: 6 mandatory columns then two alleles per variant
  PREFIX.map    classic MAP: chr, id, 0, position
  PREFIX.geno   expected converter output, one IUPAC token per sample

Encodings (`--encoding`):
  numeric   1=A 2=C 3=G 4=T, 0 = missing  (classic PED)
  letter    the bases themselves, 0 = missing

`--edge` writes a hand-built fixture instead of reading a VCF, covering the
cases the format spec calls out: both-missing, missing+observed, both bases
observed, the same base twice, a non-ACGT base, and a multi-character allele
(the converter concatenates pairs it has no IUPAC code for).
"""
import argparse
import gzip
import random
import sys

HET_PAIR = {"AC": "M", "CA": "M", "GT": "K", "TG": "K", "CT": "Y", "TC": "Y",
            "AG": "R", "GA": "R", "AT": "W", "TA": "W", "CG": "S", "GC": "S"}
NUM2BASE = {"0": "-", "1": "A", "2": "C", "3": "G", "4": "T"}
BASE2NUM = {"A": "1", "C": "2", "G": "3", "T": "4"}


def base_of(v):
    """Mirror the converter's base_of()."""
    if v is None or v == "":
        return "-"
    return NUM2BASE.get(v, v)


def merge(x, y):
    """Mirror the converter's merge_alleles()."""
    x, y = base_of(x), base_of(y)
    if x == y:
        return x
    if x == "-":
        return y
    if y == "-":
        return x
    pair = x + y
    return HET_PAIR.get(pair, pair)


def encode(base, encoding):
    """Base (or '-' for missing) -> the PED field to write."""
    if base == "-":
        return "0"
    if encoding == "numeric":
        return BASE2NUM.get(base, "0")
    return base


def gt_to_pair(gt, ref, alt):
    """VCF GT -> (allele1, allele2) as bases, '.' meaning missing."""
    a = gt.split(":")[0]
    if "/" not in a and "|" not in a:
        return "-", "-"
    sep = "|" if "|" in a else "/"
    x, y = a.split(sep)
    conv = {".": "-", "0": ref, "1": alt}
    return conv.get(x, "-"), conv.get(y, "-")


def write_vcf_case(args):
    openf = gzip.open if args.in_vcf.endswith(".gz") else open
    sites = []
    with openf(args.in_vcf, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            t = line.rstrip("\n").split("\t")
            ref, alt = t[3].upper(), t[4].upper()
            if len(ref) != 1 or len(alt) != 1:
                continue
            if args.encoding == "numeric" and (ref not in BASE2NUM or alt not in BASE2NUM):
                continue
            sites.append((t[0], t[1], ref, alt, t[9:]))
    if args.limit:
        rng = random.Random(args.seed)
        if len(sites) > args.limit:
            sites = rng.sample(sites, args.limit)
            sites.sort(key=lambda s: int(s[1]))
    return samples, sites


def write_edge_case(args):
    """Hand-built fixture: one site per merge_alleles branch."""
    samples = ["S%d" % (i + 1) for i in range(6)]
    # (chr, pos, [(a,b), ...]) with raw PED fields already in final form
    sites = [
        ("1", "100", [("1", "1"), ("1", "2"), ("2", "2"), ("0", "0"), ("0", "1"), ("4", "4")]),
        ("1", "200", [("A", "A"), ("A", "G"), ("G", "G"), ("0", "0"), ("0", "G"), ("T", "T")]),
        ("1", "300", [("0", "0"), ("0", "-"), ("-", "-"), ("N", "N"), ("N", "0"), ("0", "N")]),
        ("1", "400", [("AT", "AT"), ("AT", "G"), ("G", "AT"), ("AT", "0"), ("0", "AT"), ("AT", "-")]),
        ("1", "500", [("A", "T"), ("T", "A"), ("C", "G"), ("G", "C"), ("A", "C"), ("C", "A")]),
    ]
    return samples, sites


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--in-vcf")
    ap.add_argument("-o", "--out", required=True, help="output prefix")
    ap.add_argument("--encoding", choices=("numeric", "letter"), default="numeric")
    ap.add_argument("--edge", action="store_true", help="write the edge-case fixture")
    ap.add_argument("--limit", type=int, default=0, help="sample down to N sites")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    if args.edge:
        samples, sites = write_edge_case(args)
        rows = [(c, p, [(encode(a, args.encoding), encode(b, args.encoding))
                        for a, b in pairs]) for c, p, pairs in sites]
    else:
        if not args.in_vcf:
            ap.error("--in-vcf is required unless --edge is given")
        samples, vsites = write_vcf_case(args)
        rows = []
        for (chrom, pos, ref, alt, gts) in vsites:
            pairs = []
            for g in gts:
                a, b = gt_to_pair(g, ref, alt)
                pairs.append((encode(a, args.encoding), encode(b, args.encoding)))
            rows.append((chrom, pos, pairs))

    with open(args.out + ".ped", "w") as ped, \
         open(args.out + ".map", "w") as mp, \
         open(args.out + ".geno", "w") as geno:
        for i, s in enumerate(samples):
            ped.write("FAM%d\t%s\t0\t0\t0\t-9\t" % (i + 1, s))
            for _, _, pairs in rows:
                ped.write("%s %s " % pairs[i])
            ped.write("\n")
        for j, (chrom, pos, _) in enumerate(rows):
            mp.write("%s\tv%d\t0\t%s\n" % (chrom, j + 1, pos))
        for chrom, pos, pairs in rows:
            geno.write("%s\t%s\t%s\n" % (chrom, pos, " ".join(merge(a, b) for a, b in pairs)))

    print("wrote %s.{ped,map,geno}: %d variants x %d samples (%s)"
          % (args.out, len(rows), len(samples), args.encoding))


if __name__ == "__main__":
    sys.exit(main())
