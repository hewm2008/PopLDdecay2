#!/usr/bin/env python3
"""Generate a PLINK binary dataset (.bed/.bim/.fam) from one allele matrix, and
emit the equivalent VCF and `-InGenotype` files from that same matrix.

All three inputs describe exactly the same sites / samples / alleles, so

    PopLDdecay2 -InPlink p.bed -OutStat pl
    PopLDdecay2 -InVCF    p.vcf -OutStat vc

must produce byte-identical .stat.gz / .LD.gz / .ehh.gz (the PLINK reader
replicates the VCF site filters and 0/1/2 encoding, not the genotype ones).

PLINK .bed genotype codes, variant-major, magic 0x6c 0x1b 0x01:

    00  homozygous A1 (the first allele in .bim)
    01  missing
    10  heterozygous
    11  homozygous A2

Two modes:

  synthetic (default)  allele-frequency waves to create LD structure, plus the
                       traps that break naive implementations: sample count not
                       a multiple of 4 (byte padding), lexical chromosome order
                       (1 < 10 < 2 < X), missing calls, mono-allelic sites,
                       indels (multi-character allele), rare sites straddling
                       the -MAF threshold, and two variants sharing a position.

  --in-vcf FILE         converts a biallelic VCF (phased or unphased, gzipped or
                       not) into the PLINK triple plus a cleaned VCF of exactly
                       the sites that made it into the .bim. Records are skipped
                       when they cannot be represented identically in both
                       formats: multi-allelic ALT and non single-character
                       REF/ALT. Half-missing genotypes (0/.) are rewritten as
                       full missing calls, in the PLINK file and in the
                       cleaned VCF alike. Because the cleaned VCF is what the
                       comparison must use, the output prefix must not be the
                       input file.
"""
import argparse
import gzip
import math
import os
import random
import sys

# .bed genotype codes
HOM1, MISSING, HET, HOM2 = 0, 1, 2, 3

HET_CODE = {
    frozenset("AC"): "M", frozenset("GT"): "K", frozenset("CT"): "Y",
    frozenset("AG"): "R", frozenset("AT"): "W", frozenset("CG"): "S",
}


def iupac(a, b):
    return HET_CODE.get(frozenset((a, b)), "N")


class Variant(object):
    __slots__ = ("chrom", "pos", "vid", "a1", "a2", "gts")

    def __init__(self, chrom, pos, vid, a1, a2, gts):
        self.chrom = chrom
        self.pos = pos
        self.vid = vid
        self.a1 = a1
        self.a2 = a2
        self.gts = gts  # per-sample .bed code


# ---------------------------------------------------------------- generators

def gen_synthetic(nsamples, sites, chroms, seed):
    rng = random.Random(seed)
    bases = ["A", "C", "G", "T"]
    if chroms == 4:
        names = ["1", "2", "10", "X"]  # lexical order trap: 1 < 10 < 2 < X
    else:
        names = [str(i) for i in range(1, chroms + 1)]

    out = []
    for chrom in names:
        pos = 1000
        last = None
        for si in range(sites):
            pos += rng.randint(1, 1500)
            ref = rng.choice(bases)
            alt = rng.choice([b for b in bases if b != ref])

            if si % 50 == 0:
                p = rng.uniform(0.001, 0.004)          # rare: MAF filter territory
            else:
                p = 0.3 + 0.4 * abs(math.sin(si / 37.0))  # AF wave -> LD

            if si % 53 == 0 and si > 0:
                a1, a2 = ref + alt, ref                # indel: skipped by both paths
            elif si % 37 == 0 and si > 0:
                a1, a2 = ref, alt                      # mono-allelic below
                gts = [HOM1] * nsamples
                out.append(Variant(chrom, pos, "rs%s_%d" % (chrom, si), a1, a2, gts))
                last = out[-1]
                continue
            else:
                a1, a2 = ref, alt

            gts = []
            for i in range(nsamples):
                if out == [] and i == 0:
                    gts.append(HET)   # first data row must be a clean het:
                    continue           # the VCF path reads its phase there
                r = rng.random()
                if r < 0.04:
                    gts.append(MISSING)
                else:
                    h1 = HOM2 if rng.random() < p else HOM1
                    h2 = HOM2 if rng.random() < p else HOM1
                    gts.append(h1 if h1 == h2 else HET)
            out.append(Variant(chrom, pos, "rs%s_%d" % (chrom, si), a1, a2, gts))
            last = out[-1]

            # two variants sharing a position: the first kept occurrence wins
            if si % 71 == 0 and si > 0 and last is not None:
                out.append(Variant(chrom, pos, "rs%s_%d_dup" % (chrom, si),
                                   a2, a1, [HOM2] * (nsamples // 2) + [HET] * (nsamples - nsamples // 2)))
    return out


def open_maybe_gzip(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def convert_vcf(path):
    """VCF -> PLINK. Returns (variants, sample_names, skipped_records)."""
    variants = []
    samples = []
    skipped = 0
    with open_maybe_gzip(path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            toks = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = toks[9:]
                continue
            if len(toks) < 10:
                continue
            chrom, pos, vid, ref, alt = toks[0], toks[1], toks[2], toks[3], toks[4]
            alts = alt.split(",")
            # exactly one single-character ALT distinct from REF
            if len(alts) != 1 or len(ref) != 1 or len(alts[0]) != 1 or alts[0] == ref:
                skipped += 1
                continue
            alt = alts[0]
            gts = []
            bad = False
            for g in toks[9:]:
                a = g.split(":")[0]
                alleles = a.replace("|", "/").split("/")
                if len(alleles) != 2:
                    bad = True
                    break
                x, y = alleles
                if x == "." or y == ".":
                    # Any half-missing call becomes a full missing call in the
                    # output VCF too, so both inputs describe the same data:
                    # the VCF path reads 0/. as het-with-an-unknown-allele,
                    # which PLINK cannot express.
                    gts.append(MISSING)
                elif x == "0" and y == "0":
                    gts.append(HOM1)
                elif x == "1" and y == "1":
                    gts.append(HOM2)
                elif (x, y) in (("0", "1"), ("1", "0")):
                    gts.append(HET)
                else:
                    bad = True  # 2/3+ : multi-allelic
                    break
            if bad:
                skipped += 1
                continue
            variants.append(Variant(chrom, int(pos), vid if vid != "." else ".",
                                    ref, alt, gts))
    return variants, samples, skipped


# ------------------------------------------------------------------- writers

def write_bed(path, variants, nsamples):
    nbytes = (nsamples + 3) // 4
    with open(path, "wb") as f:
        # Three magic bytes, matching plink --make-bed: 6c 1b = variant-major,
        # 01 = SNP-only 2-bit. Writing only "6c 1b" produced a file one byte
        # shorter than any real PLINK dataset, which let the reader's own
        # off-by-one go unnoticed (both sides agreed on the wrong layout).
        f.write(b"\x6c\x1b\x01")
        for v in variants:
            buf = bytearray(nbytes)
            for k, g in enumerate(v.gts):
                buf[k >> 2] |= (g & 3) << (2 * (k & 3))
            f.write(bytes(buf))


def write_bim(path, variants):
    with open(path, "w") as f:
        for v in variants:
            f.write("%s\t%s\t0\t%d\t%s\t%s\n" % (v.chrom, v.vid, v.pos, v.a1, v.a2))


def write_fam(path, samples):
    with open(path, "w") as f:
        for i, s in enumerate(samples):
            f.write("FAM%d\t%s\t0\t0\t0\t-9\n" % (i + 1, s))


def gt_string(v):
    """VCF GT field (unphased) for one variant."""
    out = []
    for g in v.gts:
        if g == MISSING:
            out.append("./.")
        elif g == HET:
            out.append("0/1")
        elif g == HOM1:
            out.append("0/0")
        else:
            out.append("1/1")
    return "\t".join(out)


def geno_token(v, g):
    if g == MISSING:
        return "-"
    if g == HET:
        return iupac(v.a1, v.a2)
    return v.a1 if g == HOM1 else v.a2


def write_vcf(path, variants, samples):
    lines = [
        "##fileformat=VCFv4.2",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    ]
    chroms = []
    for v in variants:
        if v.chrom not in chroms:
            lines.append("##contig=<ID=%s,length=1000000000>" % v.chrom)
            chroms.append(v.chrom)
    lines.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                            "INFO", "FORMAT"] + samples))
    for v in variants:
        lines.append("\t".join([v.chrom, str(v.pos), v.vid, v.a1, v.a2, ".", "PASS",
                                ".", "GT"] + gt_string(v).split("\t")))
    data = "\n".join(lines) + "\n"
    if path.endswith(".gz"):
        with gzip.open(path, "wt") as f:
            f.write(data)
    else:
        with open(path, "w") as f:
            f.write(data)


def write_geno(path, variants, samples, with_header):
    with open(path, "w") as f:
        if with_header:
            f.write("#CHROM\tPOS\t" + " ".join(samples) + "\n")
        for v in variants:
            toks = [geno_token(v, g) for g in v.gts]
            f.write("%s\t%d\t%s\n" % (v.chrom, v.pos, " ".join(toks)))


# ---------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--out-prefix", required=True,
                    help="writes PREFIX.bed/.bim/.fam (+ PREFIX.vcf, PREFIX.geno)")
    ap.add_argument("-i", "--in-vcf", help="convert this VCF instead of generating")
    ap.add_argument("--nsamples", type=int, default=61,
                    help="samples for the synthetic mode (default 61: not a "
                         "multiple of 4, so the .bed padding bits are exercised)")
    ap.add_argument("--sites", type=int, default=1500, help="sites per chromosome")
    ap.add_argument("--chroms", type=int, default=4)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--with-header", action="store_true",
                    help="prepend the #CHROM header to the .geno file")
    args = ap.parse_args()

    p = args.out_prefix
    if args.in_vcf:
        in_vcf = os.path.abspath(args.in_vcf)
        for suffix in (".vcf", ".vcf.gz"):
            if os.path.abspath(p + suffix) == in_vcf:
                print("--out-prefix %s%s would overwrite the input VCF; "
                      "choose another prefix" % (p, suffix), file=sys.stderr)
                return 1
        variants, samples, skipped = convert_vcf(args.in_vcf)
        if not samples:
            print("no #CHROM header in %s" % args.in_vcf, file=sys.stderr)
            return 1
        if not variants:
            print("no convertible records in %s" % args.in_vcf, file=sys.stderr)
            return 1
        print("converted %s: %d sites, %d samples, %d records skipped"
              % (args.in_vcf, len(variants), len(samples), skipped))
    else:
        variants = gen_synthetic(args.nsamples, args.sites, args.chroms, args.seed)
        samples = ["S%d" % (i + 1) for i in range(args.nsamples)]

    write_bed(p + ".bed", variants, len(samples))
    write_bim(p + ".bim", variants)
    write_fam(p + ".fam", samples)
    write_vcf(p + ".vcf", variants, samples)
    write_geno(p + ".geno", variants, samples, args.with_header)
    print("wrote %s.{bed,bim,fam,vcf,geno}: %d variants, %d samples"
          % (p, len(variants), len(samples)))
    return 0


if __name__ == "__main__":
    sys.exit(main())