#!/usr/bin/env python3
"""Does a tool use the phase that is already in the VCF?

The same sites are analysed twice, once phased (0|1) and once unphased (0/1).
A tool that uses the phase answers differently for the two, a tool that
re-estimates the haplotypes from the diplotype counts answers the same.

Input (relative to the 08.Accuracy directory):
  output/poplddecay2/<TAG>.LD.gz         vs  unph_<TAG>.LD.gz
  output/plink/plink2_<TAG>.vcor          vs  unph_plink2_<TAG>.vcor

Output: output/compare/<TAG>/phase.txt

usage: phase_test.py <TAG>
"""
import gzip
import os
import sys

DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def read_pld2(path):
    out = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            out[(int(f[1]), int(f[2]))] = (float(f[5]), float(f[3]))
    return out


def read_plink2(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            k = (min(int(f[1]), int(f[3])), max(int(f[1]), int(f[3])))
            out[k] = (float(f[4]), float(f[5]))
    return out


def compare(name, a, b, lines):
    keys = sorted(set(a) & set(b))
    if not keys:
        lines.append(f"  {name:<32s} no common site pair")
        return
    same = sum(1 for k in keys if a[k] == b[k])
    dr2 = max(abs(a[k][0] - b[k][0]) for k in keys)
    ddp = max(abs(a[k][1] - b[k][1]) for k in keys)
    verdict = ("identical -> the phase is NOT used" if same == len(keys)
               else "different -> the phase IS used")
    lines.append(f"  {name:<32s} n={len(keys):<6d} identical pairs: {same:<6d} "
                 f"max|d r2|={dr2:.6f} max|d D'|={ddp:.6f}   {verdict}")


def main():
    # --full prints the banner as well; by default the screen gets the two
    # verdicts, and phase.txt keeps everything.
    full = "--full" in sys.argv
    argv = [a for a in sys.argv[1:] if not a.startswith("--")]
    tag = argv[0] if argv else "sub100_maf0.1"
    OUT = f"{DIR}/output/compare/{tag}"
    os.makedirs(OUT, exist_ok=True)
    lines = []
    lines.append("=" * 92)
    lines.append("phase test: the same sites, phased (0|1) vs unphased (0/1)")
    lines.append("=" * 92)
    tag = sys.argv[1] if len(sys.argv) > 1 else "sub100_maf0.1"
    compare("PLINK 2.0  --r2-phased",
            read_plink2(f"{DIR}/output/plink/plink2_{tag}.vcor"),
            read_plink2(f"{DIR}/output/plink/unph_plink2_{tag}.vcor"), lines)
    compare("PopLDdecay2 -OutType 8",
            read_pld2(f"{DIR}/output/poplddecay2/{tag}.LD.gz"),
            read_pld2(f"{DIR}/output/poplddecay2/unph_{tag}.LD.gz"), lines)
    lines.append("")
    lines.append("  (Haploview reads phased input only, it cannot be tested here)")
    txt = "\n".join(lines)
    with open(f"{OUT}/phase.txt", "w") as fh:
        fh.write(txt + "\n")
    if full:
        print(txt)
    else:
        for line in lines:
            if "identical pairs" in line or "Haploview reads" in line:
                print(line)
    print(f"[write] {OUT}/phase.txt")


if __name__ == "__main__":
    main()
