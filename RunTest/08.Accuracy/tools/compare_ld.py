#!/usr/bin/env python3
"""Compare the pairwise R^2 (and D', LOD, CI) of PopLDdecay2, PLINK and
Haploview against the ground truth of tools/ld_truth.py.

The headline output is the R^2 comparison:
  output/compare/sub<N>/r2_pairs.tsv    one row per site pair, the R^2 of every tool
  output/compare/sub<N>/r2_scatter.png  R^2 against each other and against the truth

Input (relative to the 08.Accuracy directory):
  output/poplddecay2/<TAG>.LD.gz      PopLDdecay2 -OutType 8
  output/plink/plink2_<TAG>.vcor       PLINK 2.0  --r2-phased
  output/haploview/<TAG>.LD            Haploview   -dprime
  output/compare/<TAG>/truth_ld.tsv    the ground truth

usage: compare_ld.py <TAG>       e.g. compare_ld.py sub100_maf0.1
"""
import gzip
import os
import statistics as st
import sys
from collections import defaultdict

import numpy as np

DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAXDIST = 300000


def out_dir(tag):
    d = f"{DIR}/output/compare/{tag}"
    os.makedirs(d, exist_ok=True)
    return d


# ------------------------------------------------------------------ readers
def read_pld2(path):
    out = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            out[(int(f[1]), int(f[2]))] = dict(
                dp=float(f[3]), lod=float(f[4]), r2=float(f[5]),
                cilow=float(f[6]), cihi=float(f[7]), dist=int(f[8]))
    return out


def read_plink19(path):
    out = {}
    with open(path) as fh:
        for i, line in enumerate(fh):
            if i == 0 or not line.strip():
                continue
            f = line.split()
            p1, p2 = int(f[1]), int(f[4])
            if abs(p2 - p1) > MAXDIST:      # the --r2 table ignores --ld-window
                continue
            out[(min(p1, p2), max(p1, p2))] = dict(
                r2=float(f[6]),
                dp=(float(f[7]) if len(f) > 7 else float("nan")),
                dist=abs(p2 - p1))
    return out


def read_plink2(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            p1, p2 = int(f[1]), int(f[3])
            if abs(p2 - p1) > MAXDIST:
                continue
            out[(min(p1, p2), max(p1, p2))] = dict(
                r2=float(f[4]), dp=float(f[5]), dist=abs(p2 - p1))
    return out


def read_haploview(path):
    out = {}
    with open(path) as fh:
        for i, line in enumerate(fh):
            if i == 0 or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            p1, p2 = int(f[0]), int(f[1])
            out[(min(p1, p2), max(p1, p2))] = dict(
                dp=float(f[2]), lod=float(f[3]), r2=float(f[4]),
                cilow=float(f[5]), cihi=float(f[6]), dist=int(f[7]))
    return out


def read_truth(path):
    out = {}
    with open(path) as fh:
        for i, line in enumerate(fh):
            if i == 0:
                continue
            f = line.rstrip("\n").split("\t")
            p1, p2 = int(f[0]), int(f[1])
            out[(min(p1, p2), max(p1, p2))] = dict(
                r2=float(f[9]), dp=float(f[8]), dist=int(f[10]),
                n=int(f[6]), n11=int(f[5]), n10=int(f[4]), n01=int(f[3]),
                n00=int(f[2]))
    return out


def read_stat(path):
    out = {}
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i == 0 or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            out[int(f[0])] = dict(mean_r2=float(f[1]), n=int(f[5]))
    return out


# ------------------------------------------------------------------ compare
def cmp_vals(a, b, label, emit, field="r2", tol=None, quiet=False):
    keys = sorted(set(a) & set(b))
    if not keys:
        emit(f"  {label:<36s} no common site pair")
        return None
    x = np.array([a[k].get(field, np.nan) for k in keys], float)
    y = np.array([b[k].get(field, np.nan) for k in keys], float)
    keep = ~np.isnan(x) & ~np.isnan(y)
    x, y = x[keep], y[keep]
    keys = [k for k, m in zip(keys, keep) if m]
    d = np.abs(x - y)
    res = dict(keys=keys, x=x, y=y, d=d,
               pear=float(np.corrcoef(x, y)[0, 1]),
               spea=float(np.corrcoef(np.argsort(np.argsort(x)),
                                      np.argsort(np.argsort(y)))[0, 1]))
    if quiet:
        return res
    extra = ""
    if tol is not None:
        extra = f"  |d|<={tol:g}: {(d <= tol).mean()*100:.2f}%"
    emit(f"  {label:<36s} n={len(keys):<6d} Pearson={res['pear']:.6f}  "
         f"Spearman={res['spea']:.6f}  mean|d|={d.mean():.6f}  "
         f"max|d|={d.max():.6f}{extra}")
    return res


def scatter(out, panels):
    """panels: list of (xlabel, ylabel, x, y, colour, title)"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, len(panels), figsize=(4.6 * len(panels), 4.4))
    if len(panels) == 1:
        axes = [axes]
    for ax, (xl, yl, x, y, col, title) in zip(axes, panels):
        ax.scatter(x, y, s=4, alpha=0.35, c=col, edgecolors="none")
        lo = float(min(np.min(x), np.min(y)))
        hi = float(max(np.max(x), np.max(y)))
        ax.plot([lo, hi], [lo, hi], "k--", lw=0.8)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        r = float(np.corrcoef(x, y)[0, 1]) if len(x) > 1 else float("nan")
        ax.set_title(f"{title}\nPearson r = {r:.6f}", fontsize=10)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)


def main():
    # --full prints the whole report to the screen as well; by default the
    # screen gets the headline table and the conclusion, and report.txt keeps
    # everything.
    full = "--full" in sys.argv
    argv = [a for a in sys.argv[1:] if not a.startswith("--")]
    tag = argv[0] if argv else "sub100_maf0.1"
    OUT = out_dir(tag)

    pld2 = read_pld2(f"{DIR}/output/poplddecay2/{tag}.LD.gz")
    p2 = read_plink2(f"{DIR}/output/plink/plink2_{tag}.vcor")
    hv = read_haploview(f"{DIR}/output/haploview/{tag}.LD")
    truth = read_truth(f"{OUT}/truth_ld.tsv")
    common = sorted(set(pld2) & set(p2) & set(hv) & set(truth))

    L = []
    A = L.append
    A("=" * 104)
    A(f"pairwise LD: PopLDdecay2(-OutType 8) vs PLINK 1.9/2.0 vs Haploview 4.2")
    A(f"data: {tag} (the first sites of base_chr22_50k.vcf.gz with MAF >= "
      "the filter of the name), 2,504 samples, window 300 kb")
    A("truth: the phased haplotypes of the VCF, counted directly (tools/ld_truth.py)")
    A("=" * 104)
    A("")
    A(f"[1] site pairs per tool (window 300 kb)")
    A(f"    PopLDdecay2 -OutType 8 : {len(pld2)}")
    A(f"    PLINK 2.0  --r2-phased  : {len(p2)}")
    A(f"    Haploview 4.2 -dprime   : {len(hv)}")
    A(f"    ground truth            : {len(truth)}")
    A(f"    common to all of them   : {len(common)}")
    A("")
    A("    sites behind those pairs. The HWE filter is switched off in all")
    A("    three (PLINK --hwe 1 1, Haploview -hwcutoff 0, PopLDdecay2 has")
    A("    none), so every program sees the same sites:")
    for nm, tbl in (("PopLDdecay2", pld2), ("PLINK 2.0 ", p2),
                    ("Haploview ", hv)):
        A(f"      {nm}: {len({x for k in tbl for x in k})}")

    A("")
    A("[2] R^2 against the ground truth  (the headline table)")
    A("    (PopLDdecay2 prints 4 decimals, Haploview 3, so they differ by rounding)")
    r_t = cmp_vals(pld2, truth, "PopLDdecay2 R^2 vs truth", A, "r2", 5e-5)
    r_h = cmp_vals(hv, truth, "Haploview  R^2 vs truth", A, "r2", 5e-4)
    r_p2 = cmp_vals(p2, truth, "PLINK 2.0  R^2 vs truth", A, "r2", 1e-3)

    A("")
    A("[3] R^2 tool against tool")
    cmp_vals(pld2, hv, "PopLDdecay2 vs Haploview", A, "r2", 5e-4)
    cmp_vals(pld2, p2, "PopLDdecay2 vs PLINK 2.0", A, "r2")
    A("")
    A("    number of site pairs where the R^2 differs by more than")
    for lab, res in (("PopLDdecay2 vs truth", r_t), ("Haploview vs truth", r_h),
                     ("PLINK 2.0 vs truth", r_p2),
                     ("PopLDdecay2 vs Haploview", cmp_vals(pld2, hv, "x", A, "r2", quiet=True)),
                     ("PopLDdecay2 vs PLINK 2.0", r_p2 and cmp_vals(pld2, p2, "x", A, "r2", quiet=True))):
        d = res["d"]
        A(f"      {lab:<28s} >0.01: {int((d > 0.01).sum()):>6d} "
          f"({(d > 0.01).mean()*100:5.2f}%)   >0.05: {int((d > 0.05).sum()):>6d} "
          f"({(d > 0.05).mean()*100:5.2f}%)   >0.20: {int((d > 0.20).sum()):>6d} "
          f"({(d > 0.20).mean()*100:5.2f}%)")

    A("")
    A("[4] D' against the ground truth")
    cmp_vals(pld2, truth, "PopLDdecay2 D' vs truth", A, "dp", 5e-5)
    cmp_vals(hv, truth, "Haploview  D' vs truth", A, "dp", 5e-4)
    cmp_vals(p2, truth, "PLINK 2.0  D' vs truth", A, "dp", 1e-3)

    A("")
    A("[5] PopLDdecay2 vs Haploview (both print D' LOD r^2 CIlow CIhi Dist)")
    cmp_vals(pld2, hv, "R^2", A, "r2", 5e-4)
    cmp_vals(pld2, hv, "D'", A, "dp", 5e-4)
    cmp_vals(pld2, hv, "LOD", A, "lod", 5e-3)
    cmp_vals(pld2, hv, "CIlow", A, "cilow", 5e-3)
    cmp_vals(pld2, hv, "CIhi", A, "cihi", 5e-3)

    A("")
    A("[6] PLINK vs truth, by how many haplotypes carry the minor allele")
    rows = []
    for k in common:
        t = truth[k]
        rows.append((t["n11"] + t["n10"] + t["n01"],
                     abs(p2[k]["r2"] - t["r2"]), abs(p2[k]["dp"] - t["dp"])))
    rows.sort()
    nb = 4
    size = max(1, len(rows) // nb)
    hdr = "mean|d D'|"
    A(f"    {'minor haplotypes':>18s} {'pairs':>7s} {'mean|d r2|':>11s} "
      f"{'max|d r2|':>10s} {hdr:>11s}")
    for b in range(nb):
        chunk = rows[b * size:(b + 1) * size] if b < nb - 1 else rows[b * size:]
        if not chunk:
            continue
        A(f"    {str(chunk[0][0])+'-'+str(chunk[-1][0]):>18s} {len(chunk):>7d} "
          f"{np.mean([c[1] for c in chunk]):>11.6f} "
          f"{np.max([c[1] for c in chunk]):>10.6f} "
          f"{np.mean([c[2] for c in chunk]):>11.6f}")

    dd = cmp_vals(pld2, p2, "x", A, "dp", quiet=True)
    A("")
    A("[7] the 10 site pairs with the largest R^2 difference "
      "(PopLDdecay2 vs PLINK 2.0)")
    A(f"    {'pos1':>10s} {'pos2':>10s} {'dist':>7s} {'PLD2 R2':>9s} "
      f"{'PLINK2 R2':>10s} {'HV R2':>9s} {'truth R2':>9s} {'PLD2 D':>8s} {'PLINK2 D':>9s}")
    for i in np.argsort(-r_p2["d"])[:10]:
        k = r_p2["keys"][i]
        A(f"    {k[0]:>10d} {k[1]:>10d} {pld2[k]['dist']:>7d} {pld2[k]['r2']:>9.4f} "
          f"{p2[k]['r2']:>10.6f} {hv[k]['r2']:>9.4f} {truth[k]['r2']:>9.6f} "
          f"{pld2[k]['dp']:>8.4f} {p2[k]['dp']:>9.4f}")

    A("")
    A("[8] LD decay curve, Mean_r^2 per distance bin (the site pairs all tools have)")
    bins = defaultdict(lambda: {"p": [], "l": [], "h": []})
    for k in common:
        b = bins[pld2[k]["dist"]]
        b["p"].append(pld2[k]["r2"])
        b["l"].append(p2[k]["r2"])
        b["h"].append(hv[k]["r2"])
    with open(f"{OUT}/curve.tsv", "w") as fh:
        fh.write("Dist\tn_pairs\tPLD2_mean_r2\tPLINK2_mean_r2\tHV_mean_r2\n")
        for b in sorted(bins):
            v = bins[b]
            fh.write(f"{b}\t{len(v['p'])}\t{st.fmean(v['p']):.6f}\t"
                     f"{st.fmean(v['l']):.6f}\t{st.fmean(v['h']):.6f}\n")
    coarse = defaultdict(lambda: {"p": [], "l": [], "h": []})
    for b, v in bins.items():
        c = coarse[(b // 10000) * 10000]
        c["p"] += v["p"]
        c["l"] += v["l"]
        c["h"] += v["h"]
    with open(f"{OUT}/curve_10kb.tsv", "w") as fh:
        fh.write("Dist_bin\tn_pairs\tPLD2_mean_r2\tPLINK2_mean_r2\tHV_mean_r2\n")
        for b in sorted(coarse):
            v = coarse[b]
            fh.write(f"{b}\t{len(v['p'])}\t{st.fmean(v['p']):.6f}\t"
                     f"{st.fmean(v['l']):.6f}\t{st.fmean(v['h']):.6f}\n")
    A(f"    {'bin':>10s} {'pairs':>7s} {'PLD2':>10s} {'PLINK2':>10s} {'Haploview':>10s}")
    for b in sorted(coarse):
        v = coarse[b]
        A(f"    {b:>10d} {len(v['p']):>7d} {st.fmean(v['p']):>10.6f} "
          f"{st.fmean(v['l']):>10.6f} {st.fmean(v['h']):>10.6f}")

    stat = read_stat(f"{DIR}/output/poplddecay2/{tag}.stat.gz")
    allp = defaultdict(list)
    for k, v in pld2.items():
        allp[v["dist"]].append(v["r2"])
    dev = [abs(st.fmean(allp[d]) - stat[d]["mean_r2"])
           for d in stat if d in allp and len(allp[d]) == stat[d]["n"]]
    if dev:
        A(f"    PopLDdecay2 .stat.gz vs its own .LD.gz, same bins: "
          f"max|diff|={max(dev):.6f} over {len(dev)} bins")

    # ---- the R^2 table -----------------------------------------------------
    with open(f"{OUT}/r2_pairs.tsv", "w") as fh:
        fh.write("pos1\tpos2\tdist\tPLD2_R2\tHV_R2\tPLINK2_R2\t"
                 "truth_R2\td_PLD2_vs_truth\td_PLD2_vs_HV\td_PLD2_vs_PLINK2\n")
        for k in common:
            p_r2 = pld2[k]["r2"]
            fh.write("\t".join(str(x) for x in [
                k[0], k[1], pld2[k]["dist"],
                f"{p_r2:.4f}", f"{hv[k]['r2']:.4f}", f"{p2[k]['r2']:.6f}",
                f"{truth[k]['r2']:.6f}",
                f"{abs(p_r2 - truth[k]['r2']):.6f}",
                f"{abs(p_r2 - hv[k]['r2']):.6f}",
                f"{abs(p_r2 - p2[k]['r2']):.6f}"]) + "\n")

    # ---- the full table ----------------------------------------------------
    with open(f"{OUT}/pairwise.tsv", "w") as fh:
        fh.write("pos1\tpos2\tdist\tPLD2_R2\tPLD2_D\tPLD2_LOD\tPLD2_CIlow\tPLD2_CIhi\t"
                 "HV_R2\tHV_D\tHV_LOD\tPLINK2_R2\tPLINK2_D\t"
                 "truth_R2\ttruth_D\tn11\tn10\tn01\tn00\n")
        for k in common:
            def g(d, f, dflt="NA"):
                return d.get(k, {}).get(f, dflt)
            fh.write("\t".join(str(x) for x in [
                k[0], k[1], pld2[k]["dist"],
                g(pld2, "r2"), g(pld2, "dp"), g(pld2, "lod"),
                g(pld2, "cilow"), g(pld2, "cihi"),
                g(hv, "r2"), g(hv, "dp"), g(hv, "lod"),
                g(p2, "r2"), g(p2, "dp"),
                g(truth, "r2"), g(truth, "dp"),
                g(truth, "n11"), g(truth, "n10"), g(truth, "n01"),
                g(truth, "n00")]) + "\n")

    # ---- the scatter plot --------------------------------------------------
    scatter(f"{OUT}/r2_scatter.png", [
        ("truth R^2", "PopLDdecay2 R^2",
         np.array([truth[k]["r2"] for k in common]),
         np.array([pld2[k]["r2"] for k in common]), "#1f77b4",
         "PopLDdecay2 vs ground truth"),
        ("PLINK 2.0 R^2", "PopLDdecay2 R^2",
         np.array([p2[k]["r2"] for k in common]),
         np.array([pld2[k]["r2"] for k in common]), "#d62728",
         "PopLDdecay2 vs PLINK 2.0"),
        ("Haploview R^2", "PopLDdecay2 R^2",
         np.array([hv[k]["r2"] for k in common]),
         np.array([pld2[k]["r2"] for k in common]), "#2ca02c",
         "PopLDdecay2 vs Haploview"),
    ])

    A("")
    A("[9] conclusion")
    A("    * PopLDdecay2 (-OutType 8) and Haploview give the same R^2, and both")
    A("      agree with the ground truth (the phased haplotypes of the VCF).")
    A("      The only differences are the printed decimals: PopLDdecay2 prints")
    A("      four, Haploview three.")
    A("    * PLINK 1.9 --r2 dprime and PLINK 2.0 --r2-phased agree with each")
    A("      other but not with the ground truth: they do not use the phase that")
    A("      is already in the VCF, they re-estimate the haplotype frequencies")
    A("      from the diplotype counts (cubic ML solution, Gaunt et al. 2007).")
    A("      The deviation is largest where the minor allele is carried by only")
    A("      a few dozen haplotypes (see [6] and [7]).")
    A("    * step4b_comparePhase.sh shows this directly: PLINK returns the same")
    A("      numbers for the phased and the unphased copy of the data.")

    txt = "\n".join(L)
    with open(f"{OUT}/report.txt", "w") as fh:
        fh.write(txt + "\n")
    if full:
        print(txt)
        for f in ("report.txt", "r2_pairs.tsv", "pairwise.tsv", "r2_scatter.png",
                  "curve.tsv", "curve_10kb.tsv", "truth_ld.tsv"):
            print(f"[write] {OUT}/{f}")
    else:
        # The headline table, and a short note instead of the whole conclusion
        # (report.txt keeps all of it). A line that starts with "[" is a
        # section header.
        brief, keep = [], False
        for line in txt.splitlines():
            if line.startswith("["):
                keep = line.startswith("[2]")
            if keep:
                brief.append(line)
        brief.append("")
        brief.append("Note:")
        brief.append("  * PopLDdecay2 and Haploview give the same r^2, and both")
        brief.append("    agree with the ground truth; only the printed decimals")
        brief.append("    differ.")
        brief.append("  * PLINK does not use the phase that is in the VCF: it")
        brief.append("    re-estimates the haplotype frequencies from the diplotype")
        brief.append("    counts. step4b_comparePhase.sh shows that directly.")
        print("\n".join(brief))
        print(f"[write] {OUT}/report.txt (and the tables next to it)")


if __name__ == "__main__":
    main()
