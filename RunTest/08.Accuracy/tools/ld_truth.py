#!/usr/bin/env python3
"""The ground truth: count the phased haplotypes of the VCF and apply the
textbook formulas for r^2 and D'.

  * every site: 1 = minor allele, 0 = major allele, 2 = missing
  * every pair: n00/n01/n10/n11 over the haplotypes that are called at both sites
  * D  = p00 - pA1 * pA2
  * r2 = D^2 / (pA1 qA1 pA2 qB2)
  * D' = |D| / Dmax

This is the same computation PopLDdecay2 does, written from scratch, so that the
three programs can be checked against something that is not one of them.

usage: ld_truth.py <in.vcf.gz> <out.tsv> [maxdist]
"""
import gzip
import sys
import numpy as np


def read_vcf(path):
    pos, codes = [], []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            pos.append(int(f[1]))
            a = []
            for g in f[9:]:
                gt = g.replace("|", "/")
                if gt[0] == ".":
                    a.extend([2, 2])
                else:
                    a.extend([int(gt[0]), int(gt[2])])
            codes.append(a)
    return np.array(pos), codes


def minor_code(arr):
    arr = np.asarray(arr)
    miss = arr == 2
    a = np.where(miss, -1, arr)
    minor = 1 if int((a == 1).sum()) < int((a == 0).sum()) else 0
    return np.where(miss, 2, (a == minor).astype(np.int8))


def pair_stats(m1, m2):
    ok = (m1 != 2) & (m2 != 2)
    n = int(ok.sum())
    if n == 0:
        return None
    a, b = m1[ok], m2[ok]
    n11 = int(((a == 1) & (b == 1)).sum())
    n10 = int(((a == 1) & (b == 0)).sum())
    n01 = int(((a == 0) & (b == 1)).sum())
    n00 = int(((a == 0) & (b == 0)).sum())
    p = np.array([n00, n01, n10, n11], dtype=float) / n
    pA1, pA2 = p[0] + p[1], p[0] + p[2]
    pB1, pB2 = 1 - pA1, 1 - pA2
    D = p[0] - pA1 * pA2
    if D >= 0:
        cal_a, cal_b = pB1 * pA2, pA1 * pB2
    else:
        D = -D
        cal_a, cal_b = pB1 * pB2, pA1 * pA2
    dmax = min(cal_a, cal_b)
    if dmax <= 0:
        return None
    return dict(n=n, n00=n00, n01=n01, n10=n10, n11=n11, D=D,
                dprime=D / dmax, r2=(D / cal_a) * (D / cal_b))


def main():
    vcf, out = sys.argv[1], sys.argv[2]
    maxdist = int(sys.argv[3]) if len(sys.argv) > 3 else 300000
    pos, codes = read_vcf(vcf)
    coded = [minor_code(c) for c in codes]
    with open(out, "w") as fh:
        fh.write("pos1\tpos2\tn00\tn01\tn10\tn11\tn_hap\tD\tdprime\tr2\tdist\n")
        n = 0
        for i in range(len(pos)):
            for j in range(i + 1, len(pos)):
                d = pos[j] - pos[i]
                if d > maxdist:
                    break
                st = pair_stats(coded[i], coded[j])
                if st is None:
                    continue
                fh.write(f"{pos[i]}\t{pos[j]}\t{st['n00']}\t{st['n01']}\t{st['n10']}"
                         f"\t{st['n11']}\t{st['n']}\t{st['D']:.6f}\t{st['dprime']:.6f}"
                         f"\t{st['r2']:.8f}\t{d}\n")
                n += 1
    print(f"[truth] {len(pos)} sites, {n} pairs -> {out}")


if __name__ == "__main__":
    main()
